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

  fn get_stroke(width:f64) -> Stroke {Stroke::new(width).with_start_cap(Cap::Butt).with_end_cap(Cap::Butt).with_join(Join::Bevel)}

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

      let dpi = 1.5;
      // Position
      let cx = 900.; let cy = 200.; let r0 = 100.;
      // Size
      let arc_len = 180; let arc_len_f = f64::from(arc_len);
      let precision_degps:f64 = 0.5; let precision_radps = precision_degps.to_radians();
      let steps_f = arc_len_f / precision_degps; //360
      let steps   = steps_f as i32;
      // Gradient / size convergence bounds
      let f_delta = 0.337; // start changing width for the first/last quarter only
      let deg_delta   	= arc_len_f *       f_delta; //18°
      let rad_delta   	= deg_delta.to_radians();
      let skip_beg_deg	= arc_len_f * (1. - f_delta);
      let skip_beg_rad	= skip_beg_deg.to_radians();
      let steps_delta 	= steps_f   * f_delta; //36
      let skip_end    	= steps_delta as i32; //36
      let skip_beg    	= steps - skip_end; //324
      // Line width
      let w1:f64 = 20.; let w2:f64 =  4.; let wavg = (w1 + w2) / 2.; // 12
      let w_delta_avg = (w2 - w1).abs() / 2.;
      let w1px = (w1 * dpi).round() / dpi; let w2px = (w2 * dpi).round() / dpi;
      let w_step = w_delta_avg / steps_delta; //12/2/45 0.13 to reach average

      // ((w1 + sign1 * w_step * r) * dpi).round() / dpi; //transition shouldn't be pixel-stepped!
      // todo3: test with dashes (unlikely to work? needs a different logic?)
        // combine vertical line with 1/2 circle
        // add dashed strokes
        // split by path
        // change stroke width in the last segments if they are covered by the 1/2 circle
          // how to calculate a match? check if we can draw a semicircle and check if overlap
        // draw previous line
        // draw splits with a different width and gradient
      let gap:f64 = 0.; // doesn't seem to 0.0001 affect anything with corrected ending style to Bevel
      let r1beg:f64 = 0.           	; let r1beg_rad = r1beg.to_radians(); //→
      let r1end = r1beg + arc_len_f	; let r1end_rad = r1end.to_radians();
      let r2beg = r1end + gap      	; let r2beg_rad = r2beg.to_radians();
      let r2end = r2beg + arc_len_f	; let r2end_rad = r2end.to_radians();
      let col_beg = css::LIME;
      let col_end = css::RED;
      let col_avg = col_beg.lerp(col_end,0.5,Default::default());

      let grad1_p0 = ( cx + r0*f64::cos( skip_beg_rad) , cy + r0*f64::sin(skip_beg_rad));
      let grad1_p1 = ( cx + r0*f64::cos( r1end_rad   ) , cy + r0*f64::sin(r1end_rad   ));
      let grad1 = Gradient::new_linear(grad1_p0, grad1_p1).with_stops([col_beg    ,col_avg]);

      let grad2_p0 = ( cx + r0*f64::cos( r2beg_rad                      ) , cy + r0*f64::sin( r2beg_rad                      ) );
      let grad2_p1 = ( cx + r0*f64::cos((r2beg + deg_delta).to_radians()) , cy + r0*f64::sin((r2beg + deg_delta).to_radians()) );
      let grad2 = Gradient::new_linear(grad2_p0, grad2_p1).with_stops([col_avg    ,col_end]);

      // Draw pre-gradwidth segment separately without the extra iterator
      let c = CircleSegment::new((cx,cy), r0,r0   ,  r1beg_rad,skip_beg_rad);
      let stroke_c = get_stroke(w1px);
      // scene.stroke(&stroke_c, Affine::IDENTITY, &col_beg, None, &c,);
      scene.stroke(&stroke_c, Affine::IDENTITY, &css::ORANGE, None, &c,); // for testing

      let steps_left_f = deg_delta / precision_degps;
      let steps_left   = steps_left_f as i32;
      let w_step_left:f64 = w_delta_avg / f64::from(steps_left); //to reach average

      let sign1 = if w1 > wavg {-1.} else if w1 < wavg {1.} else {0.}; //(to avg) ↓ if bigger, ↑ if smaller
      for i in 0..steps_left { let r = f64::from(i);
        let rad0 = (skip_beg_deg + r * precision_degps).to_radians();
        let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,precision_radps);
        //                          center  rout/in    ∠start ∠sweep
        // if i == 0 {println!("\n\n———————————————————————————————")};
        // println!("i={i}/{steps_left}  r={r:.1} r1beg={r1beg:.1} r*prec={:.1} deg={skip_beg_deg:.1} beg={:.1} end={:.1}",r * precision_degps
        //   ,       skip_beg_deg + r * precision_degps, skip_beg_deg + r * precision_degps + precision_degps);
        let cw = w1 + sign1 * r * w_step_left;
        let stroke_c = get_stroke(cw);
        scene.stroke(&stroke_c, Affine::IDENTITY, &grad1, None, &c,);
      } // ↓ in case step int conversion missed the last sliver
      let rad0_last = (skip_beg_deg + f64::from(steps_left) * precision_degps).to_radians();
      if rad0_last < r1end_rad { //println!("{rad0_last:.4}<{r1end_rad:.4} step_last {:.1}°<{:.1}° end",rad0_last.to_degrees(),r1end);
        let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0_last,r1end_rad - rad0_last);
        let stroke_c = get_stroke(w_delta_avg);
        // scene.stroke(&stroke_c, Affine::IDENTITY, &grad1, None, &c,); // use col_avg? though grad should cover
        scene.stroke(&stroke_c, Affine::IDENTITY, &css::ORANGE, None, &c,); // for testing
      }


      // Draw pos-gradwidth segment separately without the extra iterator
      let c = CircleSegment::new((cx,cy), r0,r0   ,  r2beg_rad + rad_delta, skip_beg_rad);
      // let stroke_c = get_stroke(w2px);
      let stroke_c = get_stroke(w2px).with_dashes(5.,[9.,14.]);
      // scene.stroke(&stroke_c, Affine::IDENTITY, &col_end, None, &c,);
      scene.stroke(&stroke_c, Affine::IDENTITY, &css::WHEAT, None, &c,); // for testing

      let sign2 = if w2 > wavg { 1.} else if w2 < wavg {-1.} else {0.}; //(from avg) ↑ if bigger, ↓ if smaller
      for i in 0..steps_left { let r = f64::from(i);
        let rad0 = (r2beg + r * precision_degps).to_radians();
        let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,precision_radps);
        //                          center  rout/in    ∠start ∠sweep
        // if i == 0 {println!("\n\n———————————————————————————————")};
        // println!("i={i}/{steps_left}  r={r:.1} r1beg={r1beg:.1} r*prec={:.1} deg={r2beg:.1} beg={:.1} end={:.1}",r * precision_degps
        //   ,       r2beg + r * precision_degps, r2beg + r * precision_degps + precision_degps);
        let cw = wavg + sign2 * r * w_step_left;
        let stroke_c = get_stroke(cw);
        scene.stroke(&stroke_c, Affine::IDENTITY, &grad2, None, &c,);
      } // ↓ in case step int conversion missed the last sliver
      let rad0_last = (r2beg + f64::from(steps_left) * precision_degps).to_radians();
      if rad0_last < skip_beg_rad { //println!("{rad0_last:.4} < {skip_beg_rad:.4} step_last {:.1}°<{:.1}° end",rad0_last.to_degrees(),skip_beg_deg);
        let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0_last,skip_beg_rad - rad0_last);
        let stroke_c = get_stroke(w2px);
        // scene.stroke(&stroke_c, Affine::IDENTITY, &grad2, None, &c,); // use col_avg? though grad should cover
        scene.stroke(&stroke_c, Affine::IDENTITY, &css::WHEAT, None, &c,); // for testing
      }

      // Draw debug circles showing where each gradient begins/ends
      let pstr = get_stroke(1.); // starting point bigger than the ending, angle to differentiate two curves
      let g1beg = Ellipse::new(grad1_p0, ( 2.,15.+5.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g1beg,);
      let g1end = Ellipse::new(grad1_p1, ( 2.,15.+0.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g1end,);
      let g2beg = Ellipse::new(grad2_p0, ( 2.,15.+5.),  99.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g2beg,);
      let g2end = Ellipse::new(grad2_p1, ( 2.,15.+0.),  99.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g2end,);

    }
  }
}
