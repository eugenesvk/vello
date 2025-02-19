#![cfg_attr(not(debug_assertions),allow(non_snake_case,non_upper_case_globals,non_camel_case_types))]
#![cfg_attr(    debug_assertions ,allow(non_snake_case,non_upper_case_globals,non_camel_case_types,unused_imports,unused_mut,unused_variables,dead_code,unused_assignments,unused_macros))]
use crate::{ExampleScene, SceneConfig, SceneSet};
use vello::{kurbo::{Affine, Cap},peniko::ImageQuality,};
pub fn test_scenes() -> SceneSet {test_scenes_inner()} /// All of the test scenes supported by Vello.
/// A macro which exports each passed scene indivudally This is used to avoid having to repeatedly define a
macro_rules! export_scenes {
  ($($scene_name: ident($($scene: tt)+)),*$(,)?) => {
    pub fn test_scenes_inner() -> SceneSet {
      let scenes = vec![$($scene_name()),+];
      SceneSet { scenes }
    }
    $(pub fn $scene_name() -> ExampleScene {scene!($($scene)+)})+
  };
}
/// A helper to create a shorthand name for a single scene. Used in `export_scenes`.
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
export_scenes!(stroke_styles(impls::stroke_styles(Affine::IDENTITY), "stroke_styles", false),);
/// Implementations for the test scenes. In a module because the exported [`ExampleScene`] creation functions use the same names.
mod impls {
  use std::sync::Arc;
  use crate::SceneParams;
  use kurbo::RoundedRect;
  use vello::kurbo::{Affine, BezPath, Cap, Circle, Join, PathEl, Point, Rect, Shape, Stroke, Vec2,};
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
      let join_stroke = [ MoveTo (( 0. , 0.).into()),
        CurveTo((20. , 0.).into(), (42.5, 5.).into(), ( 50., 25.).into()),
        CurveTo((57.5, 5.).into(), (80. , 0.).into(), (100.,  0.).into()),];
      let stroke1     = [ MoveTo (( 0. , 0.).into()),
        CurveTo((20. , 0.).into(), (42.5, 5.).into(), ( 50., 25.).into()),];
      let stroke2     = [ MoveTo (                   ( 50., 25.).into()),
        CurveTo((57.5, 5.).into(), (80. , 0.).into(), (100.,  0.).into()),];
      let cap_styles  = [Cap::Butt  , Cap::Round ]; //Cap::Square,
      let join_styles = [Join::Bevel, Join::Miter, Join::Round];

      // Cap and join combinations
      let mut y = 0.;
      let mut y_max: f64 = y;
      let mut color_idx = 0;
      let t = Affine::translate((60., 40.)) * Affine::scale(2.);
      y_max = y_max.max(y);
      y = 0.;
      for cap  in cap_styles  {
      for join in join_styles {
        params.text.add(scene,None,12.,None                          , Affine::translate((0., y      )) * t,
          &format!("Caps: {:?}, Joins: {:?}", cap, join),);
        scene.stroke(&Stroke::new(10.).with_caps(cap).with_join(join), Affine::translate((0., y + 30.)) * t * transform,colors[color_idx],None,
          &join_stroke);
        scene.stroke(&Stroke::new(10.).with_caps(cap).with_join(join), Affine::translate((300., y + 30.)) * t * transform,colors[color_idx],None,
          &stroke1);
        scene.stroke(&Stroke::new(20.).with_caps(cap).with_join(join), Affine::translate((300., y + 30.)) * t * transform,colors[color_idx],None,
          &stroke2);
        y += 185.;
        color_idx = (color_idx + 1) % colors.len();
      }}
    }
  }
}
