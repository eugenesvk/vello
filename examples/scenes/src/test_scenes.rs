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
  use vello::kurbo::{Affine, BezPath, Cap, Circle,CircleSegment,Ellipse, Join, PathEl, Point, Rect, Shape, Stroke, Vec2,};
  use vello::peniko::color::{palette, AlphaColor, Lch, palette::css};
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
      let dpi = 1.5;
      // Position
      let cx = 900.; let cy = 200.; let r0 = 100.;
      // Size
      let arc_len = 180; let arc_len_f = f64::from(arc_len);
      let precision_degps:f64 = 0.5; let precision_radps = precision_degps.to_radians();
      let steps_f = arc_len_f / precision_degps; //360
      let steps   = steps_f as i32;
      // Gradient / size convergence bounds
      let f_delta = 1./10.; // start changing width for the first/last quarter only
      let deg_delta  	= arc_len_f * f_delta; //18°
      let steps_delta	= steps_f   * f_delta; //36
      let skip_end   	= steps_delta as i32; //36
      let skip_beg   	= steps - skip_end; //324
      // Line width
      let w1:f64 = 20.; let w2:f64 =  4.; let wavg = (w1 + w2) / 2.; // 12
      let w1px = (w1 * dpi).round() / dpi; let w2px = (w2 * dpi).round() / dpi;
      let w_step = ((w2 - w1).abs() / 2.) / steps_delta; //12/2/45 0.13 to reach average

      let r1beg = 0. /*→*/	; let r2beg = r1beg + arc_len_f;
      let r1end = r2beg   	; let r2end = r2beg + arc_len_f;
      // todo2: check overlaps, maybe add tiny degree fractions?
      // ((w1 + sign1 * w_step * r) * dpi).round() / dpi; //transition shouldn't be pixel-stepped!
      // todo1: make the main 1st/last segment constructed in one go, no need to break into steps
      // todo2: check overlaps, maybe add tiny degree fractions?
      // todo3: test with dashes (unlikely to work? needs a different logic?)
      // todo4: add gradient
        // add same logic that gradient only starts later: so need to adjust it's point coordinates
        // make gradients also end an an average mixed color instead of stark
        //
        // NewColor = sqrt((R1^2+R2^2)/2), sqrt((G1^2+G2^2)/2), sqrt((B1^2+B2^2)/2)
        // ↑ likely bad result, use colorspace?

      let col_beg = css::LIME;
      let col_end = css::RED;
      let col_avg = col_beg.lerp(col_end,0.5,Default::default());

      let grad1_p0 = ( cx + r0*f64::cos((arc_len_f - deg_delta).to_radians()) , cy + r0*f64::sin((arc_len_f - deg_delta).to_radians()));
      let grad1_p1 = ( cx + r0*f64::cos( r1end                 .to_radians()) , cy + r0*f64::sin( r1end                   .to_radians()));
      let grad1 = Gradient::new_linear(grad1_p0, grad1_p1).with_stops([col_beg    ,col_avg]);

      let grad2_p0 = ( cx + r0*f64::cos( r2beg             .to_radians())     , cy + r0*f64::sin( r2beg             .to_radians()) );
      let grad2_p1 = ( cx + r0*f64::cos((r2beg + deg_delta).to_radians())     , cy + r0*f64::sin((r2beg + deg_delta).to_radians()) );
      // let grad2_p0 = (cx - r0                                                , cy );
      // let grad2_p1 = (cx - r0*f64::cos(             deg_delta .to_radians()) , cy - r0*f64::sin(             deg_delta .to_radians()));
      let grad2 = Gradient::new_linear(grad2_p0, grad2_p1).with_stops([col_avg    ,col_end]);

      let sign1 = if w1 > wavg {-1.} else if w1 < wavg {1.} else {0.}; //(to avg) ↓ if bigger, ↑ if smaller
      for i in 0..=steps { let r = f64::from(i); let rex = f64::from(i-skip_beg);
        let rad0 = (r1beg + r * precision_degps).to_radians();
        let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,precision_radps);
        let cw = if i > skip_beg	{w1 + sign1 * w_step * rex
        } else                  	{w1px}; //println!("wG:  {cw:.0}");
        let stroke_c = Stroke::new(cw).with_start_cap(Cap::Butt).with_end_cap(Cap::Butt);
        scene.stroke(&stroke_c, Affine::IDENTITY, &grad1, None, &c,);
      }

      let sign2 = if w2 > wavg { 1.} else if w2 < wavg {-1.} else {0.}; //(from avg) ↑ if bigger, ↓ if smaller
      for i in 0..=steps { let r = f64::from(i);
        let rad0 = (r1end + r * precision_degps).to_radians();
        let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,precision_radps);
        let cw = if i < skip_end	{wavg + sign2 * w_step * r
        } else                  	{w2px};  //println!("wR:  {cw:.0}");
        let stroke_c = Stroke::new(cw).with_start_cap(Cap::Butt).with_end_cap(Cap::Butt);
        scene.stroke(&stroke_c, Affine::IDENTITY, &grad2, None, &c,);
      }
      let pstr = Stroke::new(1.); // starting point bigger than the ending, angle to differentiate two curves
      let g1beg = Ellipse::new(grad1_p0, ( 2.,15.+5.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g1beg,);
      let g1end = Ellipse::new(grad1_p1, ( 2.,15.+0.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g1end,);
      let g2beg = Ellipse::new(grad2_p0, ( 2.,15.+5.), -66.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g2beg,);
      let g2end = Ellipse::new(grad2_p1, ( 2.,15.+0.), -66.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g2end,);


      // let cx = 900.; let cy = 200.; let r0 = 100.;
      // let grad1_p0 = (cx+(r0- 10.) , cy+ 0.);
      // let grad1_p1 = (cx-(r0-120.) , cy+140.);
      // let grad1 = Gradient::new_linear(grad1_p0, grad1_p1).with_stops([col_beg,col_end]);
      // let pstr = Stroke::new(1.);
      // let p0c = Circle ::new(grad1_p0,5.);scene.stroke(&pstr,Affine::IDENTITY, &col_beg    , None, &p0c,);
      // let p1c = Circle ::new(grad1_p1,5.);scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &p1c,);
      // let c = CircleSegment::new((cx,cy), r0,r0   ,  0.0_f64.to_radians(),180.0_f64.to_radians());
      // let stroke_c = Stroke::new(w2).with_start_cap(Cap::Butt).with_end_cap(Cap::Butt);
      // scene.stroke(&stroke_c, Affine::IDENTITY, &grad1, None, &c,);
      // // scene.stroke(&stroke_c, Affine::IDENTITY, &col_end, None, &c,);



    }
  }
}
