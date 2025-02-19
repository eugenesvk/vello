#![cfg_attr(not(debug_assertions),allow(non_snake_case,non_upper_case_globals,non_camel_case_types))]
#![cfg_attr(    debug_assertions ,allow(non_snake_case,non_upper_case_globals,non_camel_case_types,unused_imports,unused_mut,unused_variables,dead_code,unused_assignments,unused_macros))]

use crate::{ExampleScene, SceneConfig, SceneSet};
use vello::{
  kurbo::{Affine, Cap},
  peniko::ImageQuality,
};

/// All of the test scenes supported by Vello.
pub fn test_scenes() -> SceneSet {test_scenes_inner()}

/// A macro which exports each passed scene indivudally
///
/// This is used to avoid having to repeatedly define a
macro_rules! export_scenes {
  ($($scene_name: ident($($scene: tt)+)),*$(,)?) => {
    pub fn test_scenes_inner() -> SceneSet {
      let scenes = vec![
        $($scene_name()),+
      ];
      SceneSet { scenes }
    }

    $(
      pub fn $scene_name() -> ExampleScene {
        scene!($($scene)+)
      }
    )+
  };
}

/// A helper to create a shorthand name for a single scene.
/// Used in `export_scenes`.
macro_rules! scene {
  ($func:expr, $name: expr, $animated: literal) => {
    ExampleScene {
      config: SceneConfig {
        animated: $animated,
        name: $name.to_owned(),
      },
      function: Box::new($func),
    }
  };
}

export_scenes!(
  stroke_styles(impls::stroke_styles(Affine::IDENTITY), "stroke_styles", false),
);

/// Implementations for the test scenes.
/// In a module because the exported [`ExampleScene`] creation functions use the same names.
mod impls {

  use std::sync::Arc;

  use crate::SceneParams;
  use kurbo::RoundedRect;
  use vello::kurbo::{
    Affine, BezPath, Cap, Circle, Join, PathEl, Point, Rect, Shape, Stroke, Vec2,
  };
  use vello::peniko::color::{palette, AlphaColor, Lch};
  use vello::peniko::*;
  use vello::*;

  pub(super) fn stroke_styles(transform: Affine) -> impl FnMut(&mut Scene, &mut SceneParams<'_>) {
    use PathEl::*;
    move |scene, params| {
      let colors = [
        Color::from_rgb8(140, 181, 236),
        Color::from_rgb8(246, 236, 202),
        Color::from_rgb8(201, 147, 206),
        Color::from_rgb8(150, 195, 160),
      ];
      let simple_stroke = [MoveTo((0., 0.).into()), LineTo((100., 0.).into())];
      let join_stroke = [
        MoveTo (( 0. , 0.).into()),
        CurveTo((20. , 0.).into(), (42.5, 5.).into(), ( 50., 25.).into()),
        CurveTo((57.5, 5.).into(), (80. , 0.).into(), (100.,  0.).into()),
      ];
      let miter_stroke = [
        MoveTo((0., 0.).into()),
        LineTo((90., 16.).into()),
        LineTo((0., 31.).into()),
        LineTo((90., 46.).into()),
      ];
      let closed_strokes = [
        MoveTo((0., 0.).into()),
        LineTo((90., 21.).into()),
        LineTo((0., 42.).into()),
        ClosePath,
        MoveTo((200., 0.).into()),
        CurveTo((100., 72.).into(), (300., 72.).into(), (200., 0.).into()),
        ClosePath,
        MoveTo((290., 0.).into()),
        CurveTo((200., 72.).into(), (400., 72.).into(), (310., 0.).into()),
        ClosePath,
      ];
      let cap_styles = [Cap::Butt, Cap::Square, Cap::Round];
      let join_styles = [Join::Bevel, Join::Miter, Join::Round];
      let miter_limits = [4., 6., 0.1, 10.];

      // Simple strokes with cap combinations
      let t = Affine::translate((60., 40.)) * Affine::scale(2.);
      let mut y = 0.;
      let mut color_idx = 0;
      for start in cap_styles {
        for end in cap_styles {
          params.text.add(scene,None,12.,None,
            Affine::translate((0., y)) * t,
            &format!("Start cap: {:?}, End cap: {:?}", start, end),
          );
          scene.stroke(
            &Stroke::new(20.).with_start_cap(start).with_end_cap(end),
            Affine::translate((0., y + 30.)) * t * transform,
            colors[color_idx],
            None,
            &simple_stroke,
          );
          y += 180.;
          color_idx = (color_idx + 1) % colors.len();
        }
      }
      // Dashed strokes with cap combinations
      let t = Affine::translate((450., 0.)) * t;
      let mut y_max = y;
      y = 0.;
      for start in cap_styles {
        for end in cap_styles {
          params.text.add(scene,None,12.,None,
            Affine::translate((0., y)) * t,
            &format!("Dashing - Start cap: {:?}, End cap: {:?}", start, end),
          );
          scene.stroke(
            &Stroke::new(20.)
              .with_start_cap(start)
              .with_end_cap(end)
              .with_dashes(0., [10.0, 21.0]),
            Affine::translate((0., y + 30.)) * t * transform,
            colors[color_idx],
            None,
            &simple_stroke,
          );
          y += 180.;
          color_idx = (color_idx + 1) % colors.len();
        }
      }

      // Cap and join combinations
      let t = Affine::translate((550., 0.)) * t;
      y_max = y_max.max(y);
      y = 0.;
      for cap in cap_styles {
        for join in join_styles {
          params.text.add(scene,None,12.,None,
            Affine::translate((0., y)) * t,
            &format!("Caps: {:?}, Joins: {:?}", cap, join),
          );
          scene.stroke(
            &Stroke::new(20.).with_caps(cap).with_join(join),
            Affine::translate((0., y + 30.)) * t * transform,
            colors[color_idx],
            None,
            &join_stroke,
          );
          y += 185.;
          color_idx = (color_idx + 1) % colors.len();
        }
      }

      // Miter limit
      let t = Affine::translate((500., 0.)) * t;
      y_max = y_max.max(y);
      y = 0.;
      for ml in miter_limits {
        params.text.add(scene,None,12.,None,
          Affine::translate((0., y)) * t,
          &format!("Miter limit: {}", ml),
        );
        scene.stroke(
          &Stroke::new(10.)
            .with_caps(Cap::Butt)
            .with_join(Join::Miter)
            .with_miter_limit(ml),
          Affine::translate((0., y + 30.)) * t * transform,
          colors[color_idx],
          None,
          &miter_stroke,
        );
        y += 180.;
        color_idx = (color_idx + 1) % colors.len();
      }

      // Closed paths
      for (i, join) in join_styles.iter().enumerate() {
        params.text.add(scene,None,12.,None,
          Affine::translate((0., y)) * t,
          &format!("Closed path with join: {:?}", join),
        );
        // The cap style is not important since a closed path shouldn't have any caps.
        scene.stroke(
          &Stroke::new(10.)
            .with_caps(cap_styles[i])
            .with_join(*join)
            .with_miter_limit(5.),
          Affine::translate((0., y + 30.)) * t * transform,
          colors[color_idx],
          None,
          &closed_strokes,
        );
        y += 180.;
        color_idx = (color_idx + 1) % colors.len();
      }
      y_max = y_max.max(y);
      // The closed_strokes has a maximum x of 400, `t` has a scale of `2.`
      let x_max = t.translation().x + 400. * 2. + 50.; // Give 50px of padding to account for `transform`
      params.resolution = Some((x_max, y_max).into());
    }
  }
}
