// Copyright 2022 the Vello Authors
// SPDX-License-Identifier: Apache-2.0 OR MIT

//! Winit example.

// The following lints are part of the Linebender standard set,
// but resolving them has been deferred for now.
// Feel free to send a PR that solves one or more of these.
#![allow(
  unreachable_pub,
  clippy::allow_attributes_without_reason,
  clippy::cast_possible_truncation,
  clippy::shadow_unrelated
)]

use std::num::NonZeroUsize;
use std::sync::Arc;

use std::time::Instant;
use vello::low_level::DebugLayers;
use winit::application::ApplicationHandler;
use winit::event::{ElementState, MouseButton, MouseScrollDelta, TouchPhase, WindowEvent};
use winit::keyboard::{Key, ModifiersState, NamedKey};


use clap::Parser;
use scenes::{ExampleScene, ImageCache, SceneParams, SceneSet, SimpleText};
use vello::kurbo::{Affine, Vec2};
use vello::peniko::{color::palette, Color};
use vello::util::{RenderContext, RenderSurface};
use vello::{low_level::BumpAllocators, AaConfig, Renderer, RendererOptions, Scene};

use winit::dpi::LogicalSize;
use winit::event_loop::EventLoop;
use winit::window::{Window, WindowAttributes};

use vello::wgpu;

mod hot_reload;

#[derive(Parser, Debug)]
#[command(about, long_about = None, bin_name="cargo run -p with_winit --")]
struct Args {
  /// Which scene (index) to start on Switch between scenes with left and right arrow keys
  #[arg(long)]
  scene: Option<i32>,
  #[command(flatten)]
  args: scenes::Arguments,
  #[arg(long)]
  /// Whether to use CPU shaders
  use_cpu: bool,
  /// Used to disable vsync at startup. Can be toggled with the "V" key.
  ///
  /// This setting is useful for Android, where it might be harder to press this key
  #[arg(long)]
  startup_vsync_off: bool,
  /// Used to enable gpu profiling at startup. Can be toggled with the "G" key
  ///
  /// It is off by default because it has adverse performance characteristics
  #[arg(long)]
  #[cfg(feature = "wgpu-profiler")]
  startup_gpu_profiling_on: bool,
  /// Whether to force initialising the shaders serially (rather than spawning threads)
  /// This has no effect on wasm, and defaults to 1 on macOS for performance reasons
  ///
  /// Use `0` for an automatic choice
  #[arg(long, default_value_t=default_threads())]
  num_init_threads: usize,
  /// Use the asynchronous pipeline (if available) for rendering
  ///
  /// The asynchronous pipeline is one approach for robust memory - see
  /// <https://github.com/linebender/vello/issues/366>
  ///
  /// However, it also has potential latency issues, especially for
  /// accessibility technology, as it (currently) blocks the main thread for
  /// extended periods
  #[arg(long)]
  async_pipeline: bool,
}

fn default_threads() -> usize {
  #[cfg    (target_os = "macos") ]return 1;
  #[cfg(not(target_os = "macos"))]return 0;
}

struct RenderState<'s> {
  // SAFETY: We MUST drop the surface before the `window`, so the fields
  // must be in this order
  surface: RenderSurface<'s>,
  window: Arc<Window>,
}

// TODO: Make this set configurable through the command line
// Alternatively, load anti-aliasing shaders on demand/asynchronously
const AA_CONFIGS: [AaConfig; 3] = [AaConfig::Area, AaConfig::Msaa8, AaConfig::Msaa16];

struct VelloApp<'s> {
  context: RenderContext,
  renderers: Vec<Option<Renderer>>,
  state: Option<RenderState<'s>>,
  // Whilst suspended, we drop `render_state`, but need to keep the same window.
  // If render_state exists, we must store the window in it, to maintain drop order
  cached_window: Option<Arc<Window>>,

  use_cpu: bool,
  num_init_threads: usize,

  scenes: Vec<ExampleScene>,
  scene: Scene,
  fragment: Scene,
  simple_text: SimpleText,
  images: ImageCache,

  base_color: Option<Color>,
  async_pipeline: bool,

  scene_complexity: Option<BumpAllocators>,

  complexity_shown: bool,
  vsync_on: bool,

  #[cfg(feature = "wgpu-profiler")]
  gpu_profiling_on: bool,
  #[cfg(feature = "wgpu-profiler")]
  profile_stored: Option<Vec<wgpu_profiler::GpuTimerQueryResult>>,

  // We allow cycling through AA configs in either direction, so use a signed index
  aa_config_ix: i32,

  frame_start_time: Instant,
  start: Instant,

  transform: Affine,
  mouse_down: bool,
  prior_position: Option<Vec2>,
  scene_ix: i32, // We allow looping left and right through the scenes, so use a signed index
  complexity: usize,

  prev_scene_ix: i32,
  modifiers: ModifiersState,

  debug: DebugLayers,
}

impl ApplicationHandler<UserEvent> for VelloApp<'_> {
  fn resumed(&mut self, event_loop: &winit::event_loop::ActiveEventLoop) {
    let None = self.state else {
      return;
    };
    let window = self
      .cached_window
      .take()
      .unwrap_or_else(|| Arc::new(event_loop.create_window(window_attributes()).unwrap()));
    let size = window.inner_size();
    let present_mode = if self.vsync_on {wgpu::PresentMode::AutoVsync
    } else {wgpu::PresentMode::AutoNoVsync};
    let surface_future =
      self.context.create_surface(window.clone(), size.width, size.height, present_mode);
    // We need to block here, in case a Suspended event appeared
    let surface = pollster::block_on(surface_future).expect("Error creating surface");
    self.state = {
      let render_state = RenderState { window, surface };
      self.renderers.resize_with(self.context.devices.len(), || None);
      let id = render_state.surface.dev_id;
      self.renderers[id].get_or_insert_with(|| {
        let start = Instant::now();
        let renderer = Renderer::new(
          &self.context.devices[id].device,
          RendererOptions {
            use_cpu: self.use_cpu,
            antialiasing_support: AA_CONFIGS.iter().copied().collect(),
            num_init_threads: NonZeroUsize::new(self.num_init_threads),
          },
        ).map_err(|e| {anyhow::format_err!("{e}")}) // Pretty-print any renderer creation error using Display formatting before unwrapping.
        .expect("Failed to create renderer");
        log::info!("Creating renderer {id} took {:?}", start.elapsed());
        #[cfg(feature = "wgpu-profiler")]
        let mut renderer = renderer;
        #[cfg(feature = "wgpu-profiler")]
        renderer
          .profiler
          .change_settings(wgpu_profiler::GpuProfilerSettings {
            enable_timer_queries: self.gpu_profiling_on,
            enable_debug_groups: self.gpu_profiling_on,
            ..Default::default()
          })
          .expect("Not setting max_num_pending_frames");
        renderer
      });
      Some(render_state)
    };
  }

  fn window_event(
    &mut self,
    event_loop: &winit::event_loop::ActiveEventLoop,
    window_id: winit::window::WindowId,
    event: WindowEvent,
  ) {
    let Some(render_state) = &mut self.state else {
      return;
    };
    if render_state.window.id() != window_id {
      return;
    }
    let _span = if !matches!(event, WindowEvent::RedrawRequested) {
      Some(tracing::trace_span!("Handling window event", ?event).entered())
    } else {
      None
    };
    match event {
      WindowEvent::CloseRequested => event_loop.exit(),
      WindowEvent::ModifiersChanged(m) => self.modifiers = m.state(),
      WindowEvent::KeyboardInput { event, .. } => {
        if event.state == ElementState::Pressed {
          match event.logical_key.as_ref() {
            Key::Named(NamedKey::ArrowLeft) => {
              println!("← {}",self.scene_ix);
              self.scene_ix = self.scene_ix.saturating_sub(1);
            }
            Key::Named(NamedKey::ArrowRight) => {
              println!("→ {}",self.scene_ix);
              self.scene_ix = self.scene_ix.saturating_add(1);
              println!("→ {}",self.scene_ix);
            }
            Key::Named(NamedKey::ArrowUp) => self.complexity += 1,
            Key::Named(NamedKey::ArrowDown) => {
              self.complexity = self.complexity.saturating_sub(1);
            }
            Key::Named(NamedKey::Space) => {
              self.transform = Affine::IDENTITY;
            }
            Key::Character(char) => {
              // TODO: Have a more principled way of handling modifiers on keypress
              // see e.g. https://xi.zulipchat.com/#narrow/stream/351333-glazier/topic/Keyboard.20shortcuts
              let char = char.to_lowercase();
              match char.as_str() {
                "q" | "e" => {
                  if let Some(prior_position) = self.prior_position {
                    let is_clockwise = char == "e";
                    let angle = if is_clockwise { -0.05 } else { 0.05 };
                    self.transform = Affine::translate(prior_position)
                      * Affine::rotate(angle)
                      * Affine::translate(-prior_position)
                      * self.transform;
                  }
                }
                "d" => {self.complexity_shown = !self.complexity_shown;}
                "m" => {
                  self.aa_config_ix = if self.modifiers.shift_key() {self.aa_config_ix.saturating_sub(1)
                  } else {self.aa_config_ix.saturating_add(1)
                  };}
                #[cfg(feature = "wgpu-profiler")]
                "p" => {if let Some(renderer) =&self.renderers[render_state.surface.dev_id]{
                  store_profiling(renderer, &self.profile_stored);}
                }
                #[cfg(feature = "wgpu-profiler")]
                "g" => {
                  self.gpu_profiling_on = !self.gpu_profiling_on;
                  if let Some(renderer) = &mut self.renderers[render_state.surface.dev_id] {
                    renderer.profiler.change_settings(wgpu_profiler::GpuProfilerSettings {
                        enable_timer_queries: self.gpu_profiling_on,
                        enable_debug_groups : self.gpu_profiling_on,
                        ..Default::default()
                      })
                      .expect("Not setting max_num_pending_frames");
                  }
                }
                "v" => {
                  self.vsync_on = !self.vsync_on;
                  self.context.set_present_mode(
                    &mut render_state.surface,
                    if self.vsync_on {
                      wgpu::PresentMode::AutoVsync
                    } else {
                      wgpu::PresentMode::AutoNoVsync
                    },
                  );
                }
                debug_layer @ ("1" | "2" | "3" | "4") => {
                  match debug_layer {
                    "1" => {
                      self.debug.toggle(DebugLayers::BOUNDING_BOXES);
                    }
                    "2" => {
                      self.debug.toggle(DebugLayers::LINESOUP_SEGMENTS);
                    }
                    "3" => {
                      self.debug.toggle(DebugLayers::LINESOUP_POINTS);
                    }
                    "4" => {
                      self.debug.toggle(DebugLayers::VALIDATION);
                    }
                    _ => unreachable!(),
                  }
                  if !self.debug.is_empty() && !self.async_pipeline {
                    log::warn!("Debug Layers won't work without using `--async-pipeline`. Requested {:?}", self.debug);
                  }
                }
                _ => {}
              }
            }
            Key::Named(NamedKey::Escape) => event_loop.exit(),
            _ => {}
          }
        }
      }
      WindowEvent::Touch(touch) => {
        match touch.phase {
          TouchPhase::Started => {}
          TouchPhase::Ended | TouchPhase::Cancelled => {}
          TouchPhase::Moved => (),
        }
      }
      WindowEvent::Resized(size) => {
        if let Some(RenderState { surface, window }) = &mut self.state {
          self.context.resize_surface(surface, size.width, size.height);
          window.request_redraw();
        }
      }
      WindowEvent::MouseInput { state, button, .. } => {
        if button == MouseButton::Left {
          self.mouse_down = state == ElementState::Pressed;
        }
      }
      WindowEvent::MouseWheel { delta, .. } => {
        const BASE: f64 = 1.05;
        const PIXELS_PER_LINE: f64 = 20.0;

        if let Some(prior_position) = self.prior_position {
          let exponent = if let MouseScrollDelta::PixelDelta(delta) = delta {
            delta.y / PIXELS_PER_LINE
          } else if let MouseScrollDelta::LineDelta(_, y) = delta {
            y as f64
          } else {
            0.0
          };
          self.transform = Affine::translate(prior_position)
            * Affine::scale(BASE.powf(exponent))
            * Affine::translate(-prior_position)
            * self.transform;
        } else {
          log::warn!("Scrolling without mouse in window; this shouldn't be possible");
        }
      }
      WindowEvent::CursorLeft { .. } => {
        self.prior_position = None;
      }
      WindowEvent::CursorMoved { position, .. } => {
        let position = Vec2::new(position.x, position.y);
        if self.mouse_down {
          if let Some(prior) = self.prior_position {
            self.transform = Affine::translate(position - prior) * self.transform;
          }
        }
        self.prior_position = Some(position);
      }
      WindowEvent::RedrawRequested => {
        let _rendering_span = tracing::trace_span!("Actioning Requested Redraw").entered();
        let encoding_span = tracing::trace_span!("Encoding scene").entered();

        render_state.window.request_redraw();

        let Some(RenderState { surface, window }) = &self.state else {
          return;
        };
        let width = surface.config.width;
        let height = surface.config.height;
        let device_handle = &self.context.devices[surface.dev_id];

        // Allow looping forever
        self.scene_ix = self.scene_ix.rem_euclid(self.scenes.len() as i32);
        self.aa_config_ix = self.aa_config_ix.rem_euclid(AA_CONFIGS.len() as i32);

        let example_scene = &mut self.scenes[self.scene_ix as usize];
        if self.prev_scene_ix != self.scene_ix {
          self.transform = Affine::IDENTITY;
          self.prev_scene_ix = self.scene_ix;
          window.set_title(&format!("Vello demo - {}", example_scene.config.name));
        }
        self.fragment.reset();

        let mut scene_params = SceneParams { // !!!
          time: self.start.elapsed().as_secs_f64(),
          text: &mut self.simple_text,
          images: &mut self.images,
          resolution: None,
          base_color: None,
          interactive: true,
          complexity: self.complexity,
        };
        example_scene.function.render(&mut self.fragment, &mut scene_params);

        // If the user specifies a base color in the CLI we use that. Otherwise we use any color specified by the scene. The default is black.
        let base_color = self.base_color.or(scene_params.base_color).unwrap_or(palette::css::BLACK);
        let antialiasing_method = AA_CONFIGS[self.aa_config_ix as usize];
        let render_params = vello::RenderParams {
          base_color,
          width,
          height,
          antialiasing_method,
        };
        self.scene.reset();
        let mut transform = self.transform;
        if let Some(resolution) = scene_params.resolution {
          // Automatically scale the rendering to fill as much of the window as possible TODO: Apply svg view_box, somehow
          let factor = Vec2::new(width as f64, height as f64);
          let scale_factor = (factor.x / resolution.x).min(factor.y / resolution.y);
          transform *= Affine::scale(scale_factor);
        }

        self.scene.append(&self.fragment, Some(transform));
        drop(encoding_span);
        let render_span = tracing::trace_span!("Dispatching render").entered();
        // Note: we don't run the async/"robust" pipeline on web, as
        // it requires more async wiring for the readback. See
        // [#gpu > async on wasm](https://xi.zulipchat.com/#narrow/stream/197075-gpu/topic/async.20on.20wasm)
        #[expect(
          deprecated,
          reason = "We still want to use the async pipeline for the debug layers"
        )]
        if self.async_pipeline {
          self.scene_complexity = vello::util::block_on_wgpu(
            &device_handle.device,
            self.renderers[surface.dev_id]
              .as_mut()
              .unwrap()
              .render_to_texture_async(
                &device_handle.device,
                &device_handle.queue,
                &self.scene,
                &surface.target_view,
                &render_params,
                self.debug,
              ),
          ).expect("failed to render to texture");
        } else {
          self.renderers[surface.dev_id].as_mut().unwrap().render_to_texture(
            &device_handle.device,
            &device_handle.queue,
            &self.scene,
            &surface.target_view,
            &render_params,
          ).expect("failed to render to texture");
        }
        drop(render_span);

        let texture_span = tracing::trace_span!("Blitting to surface").entered();
        let surface_texture = surface
          .surface
          .get_current_texture()
          .expect("failed to get surface texture");
        // Perform the copy
        // (TODO: Does it improve throughput to acquire the surface after the previous texture render has happened?)
        let mut encoder =
          device_handle
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
              label: Some("Surface Blit"),
            });
        surface.blitter.copy(
          &device_handle.device,
          &mut encoder,
          &surface.target_view,
          &surface_texture
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default()),
        );
        device_handle.queue.submit([encoder.finish()]);
        surface_texture.present();
        drop(texture_span);

        {
          let _poll_span = tracing::trace_span!("Polling wgpu device").entered();
          device_handle.device.poll(wgpu::Maintain::Poll);
        }
        let new_time = Instant::now();
        self.frame_start_time = new_time;
      }
      _ => {}
    }
  }

  fn about_to_wait(&mut self, _event_loop: &winit::event_loop::ActiveEventLoop) {
    if let Some(render_state) = &mut self.state {
      render_state.window.request_redraw();
    }
  }

  fn user_event(&mut self, _event_loop: &winit::event_loop::ActiveEventLoop, event: UserEvent) {
    match event {
      UserEvent::HotReload => {
        let Some(render_state) = &mut self.state else {
          return;
        };
        let device_handle = &self.context.devices[render_state.surface.dev_id];
        log::info!("==============\nReloading shaders");
        let start = Instant::now();
        let result = self.renderers[render_state.surface.dev_id]
          .as_mut()
          .unwrap()
          .reload_shaders(&device_handle.device);
        // We know that the only async here (`pop_error_scope`) is actually sync, so blocking is fine
        match pollster::block_on(result) {
          Ok(_) => log::info!("Reloading took {:?}", start.elapsed()),
          Err(e) => log::error!("Failed to reload shaders: {e}"),
        }
      }
    }
  }

  fn suspended(&mut self, _event_loop: &winit::event_loop::ActiveEventLoop) {
    log::info!("Suspending");
    // When we suspend, we need to remove the `wgpu` Surface
    if let Some(render_state) = self.state.take() {
      self.cached_window = Some(render_state.window);
    }
  }
}

fn run(
  event_loop: EventLoop<UserEvent>,
  args: Args,
  scenes: SceneSet,
  render_cx: RenderContext,
) {
  use winit::keyboard::ModifiersState;

  let (render_state, renderers) = (None::<RenderState<'_>>, vec![]);

  // The design of `RenderContext` forces delayed renderer initialisation to
  // not work on wasm, as WASM futures effectively must be 'static.
  // Otherwise, this could work by sending the result to event_loop.proxy
  // instead of blocking

  let debug = DebugLayers::none();

  let mut app = VelloApp {
    context: render_cx,
    renderers,
    state: render_state,
    cached_window: None,
    use_cpu: args.use_cpu,
    num_init_threads: args.num_init_threads,
    scenes: scenes.scenes,
    scene: Scene::new(),
    fragment: Scene::new(),
    simple_text: SimpleText::new(),
    images: ImageCache::new(),
    base_color: args.args.base_color,
    async_pipeline: args.async_pipeline,
    scene_complexity: None,
    complexity_shown: false,
    vsync_on: !args.startup_vsync_off,

    #[cfg(feature = "wgpu-profiler")]
    gpu_profiling_on: args.startup_gpu_profiling_on,
    #[cfg(feature = "wgpu-profiler")]
    profile_stored: None,

    aa_config_ix: 0,

    frame_start_time: Instant::now(),
    start: Instant::now(),

    transform: Affine::IDENTITY,
    mouse_down: false,
    prior_position: None,
    scene_ix: args.scene.unwrap_or(0),
    complexity: 0,
    prev_scene_ix: 0,
    modifiers: ModifiersState::default(),
    debug,
  };

  event_loop.run_app(&mut app).expect("run to completion");
}

#[cfg(feature = "wgpu-profiler")]
/// A function extracted to fix rustfmt
fn store_profiling(
  renderer: &Renderer,
  profile_stored: &Option<Vec<wgpu_profiler::GpuTimerQueryResult>>,
) {
  if let Some(profile_result) = &renderer.profile_result.as_ref().or(profile_stored.as_ref()) {
    if !profile_result.is_empty() { // There can be empty results if the required features aren't supported
      let path = std::path::Path::new("trace.json");
      match wgpu_profiler::chrometrace::write_chrometrace(path, profile_result) {
        Ok(()) => {println!("Wrote trace to path {path:?}");}
        Err(e) => {log::warn!("Failed to write trace {e}");}
      }
    }
  }
}

fn window_attributes() -> WindowAttributes {
  Window::default_attributes().with_inner_size(LogicalSize::new(800,620)).with_resizable(true).with_title("Vello demo")
}
#[derive(Debug)]
enum UserEvent {HotReload,}

/// Entry point.
pub fn main() -> anyhow::Result<()> {
  // TODO: initializing both env_logger and console_logger fails on wasm.
  // Figure out a more principled approach.
  env_logger::builder()
    .format_timestamp(Some(env_logger::TimestampPrecision::Millis))
    .filter_level(log::LevelFilter::Warn)
    .init();
  let args = parse_arguments();
  let scenes = args.args.select_scene_set()?;
  if let Some(scenes) = scenes {
    println!("scenes");
    let event_loop = EventLoop::<UserEvent>::with_user_event().build()?;
    let render_cx = RenderContext::new();
    let proxy = event_loop.create_proxy();
    let _keep = hot_reload::hot_reload(move || {proxy.send_event(UserEvent::HotReload).ok().map(drop)});
    run(event_loop, args, scenes, render_cx);
  }
  Ok(())
}

fn parse_arguments() -> Args {
  // We allow baking in arguments at compile time. This is especially useful for Android and WASM.
  // This is used on desktop platforms to allow debugging the same settings
  let args = if let Some(args) = option_env!("VELLO_STATIC_ARGS") {
    // We split by whitespace here to allow passing multiple arguments
    // In theory, we could do more advanced parsing/splitting (e.g. using quotes),
    // but that would require a lot more effort

    // We `chain` in a fake binary name, because clap ignores the first argument otherwise
    // Ideally, we'd use the `no_binary_name` argument, but setting that at runtime would
    // require globals or some worse hacks
    Args::parse_from(std::iter::once("with_winit").chain(args.split_ascii_whitespace()))
  } else {
    Args::parse()
  };
  args
}

#[cfg(all(feature = "_ci_dep_features_to_test", test))]
#[test]
// This just tests that the "kurbo" dependency we enable schemars for
// aligns to the same version that vello's peniko dependency resolves to.
fn test_kurbo_schemars_with_peniko() {
  use std::marker::PhantomData;
  #[expect(unused_qualifications)]
  let _: PhantomData<kurbo::Rect> = PhantomData::<vello::peniko::kurbo::Rect>;
}
