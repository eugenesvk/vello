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
  use vello::kurbo::{Affine, BezPath, Cap, Circle,CircleSegment,Ellipse,Line,Join, PathEl, Point, Rect, Shape, Stroke, Vec2,};
  use vello::peniko::color::{palette, AlphaColor, Lch, palette::css};
  use vello::peniko::*;
  use vello::*;
  use std::f64::consts as f64c;

  fn get_stroke    (width:f64) -> Stroke {Stroke::new(width).with_start_cap(Cap::Butt).with_end_cap(Cap::Round).with_join(Join::Bevel)}
  fn get_stroke_end(width:f64) -> Stroke {Stroke::new(width).with_start_cap(Cap::Butt).with_end_cap(Cap::Butt ).with_join(Join::Bevel)}
  fn get_stroke1   (width:f64) -> Stroke {Stroke::new(width).with_start_cap(Cap::Butt).with_end_cap(Cap::Square).with_join(Join::Round)}

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
      let cx = 900.; let cy = 200.; let r0 = 95.5; //600 circum len 300 half
      // Size
      let arc_len = 180; let arc_len_f = f64::from(arc_len);
      let precision_degps:f64 = 0.5; let precision_radps = precision_degps.to_radians(); //0.00873
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

      // TODO: change all circle segments to .outer_arc() otherwise we're drawing 2! lines
        // but this introduces strange artifacts :(( joins are visible, changing ending to round fixes it, but i want sharp ends!
        // ask how to fix them?
      // ((w1 + sign1 * w_step * r) * dpi).round() / dpi; //transition shouldn't be pixel-stepped!
      // todo3: test with dashes (unlikely to work? needs a different logic?)
        // Approach1:
          // calculate dashed sequence [1,2,1,4] length in radians (1+2+1+4=8px == 1 rad at radius X)
          // for each iterated gradwidht tiny segment, get remainder of divis by this iter length
          // then determine, which sublength this short line belongs to (2nd: from 1/8 to (1+3)/8)
          // then if this is an odd "dash" segment, draw, if not, skip
        // Approach2: failed
          // combine vertical line with 1/2 circle
          // add dashed strokes
          // split by path
          // change stroke width in the last segments if they are covered by the 1/2 circle
            // not possible???? how to calculate a match? check if we can draw a semicircle and check if overlap
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
      let c = CircleSegment::new((cx,cy), r0,0.   ,  r1beg_rad,skip_beg_rad).outer_arc();
      let stroke_c = get_stroke_end(w1px);
      // scene.stroke(&stroke_c, Affine::IDENTITY, &col_beg, None, &c,);
      scene.stroke(&stroke_c, Affine::IDENTITY, &css::ORANGE, None, &c,); // for testing

      let steps_left_f = deg_delta / precision_degps;
      let steps_left   = steps_left_f as i32;
      let w_step_left:f64 = w_delta_avg / f64::from(steps_left); //to reach average

      let sign1 = if w1 > wavg {-1.} else if w1 < wavg {1.} else {0.}; //(to avg) ↓ if bigger, ↑ if smaller
      for i in 0..steps_left { let r = f64::from(i);
        let rad0 = (skip_beg_deg + r * precision_degps).to_radians();
        let c = CircleSegment::new((cx,cy), r0,0.   ,  rad0,precision_radps).outer_arc();
        //                          center  rout/in    ∠start ∠sweep
        // if i == 0 {println!("\n\n———————————————————————————————")};
        // println!("i={i}/{steps_left}  r={r:.1} r1beg={r1beg:.1} r*prec={:.1} deg={skip_beg_deg:.1} beg={:.1} end={:.1}",r * precision_degps
        //   ,       skip_beg_deg + r * precision_degps, skip_beg_deg + r * precision_degps + precision_degps);
        let cw = w1 + sign1 * r * w_step_left;
        let stroke_c = if i >=  (steps_left - 2) {get_stroke_end(cw) //last steps no round ends
        } else {get_stroke(cw)}; // round ends to fix the outer arc artifacts when joining in rectangles without curves
        scene.stroke(&stroke_c, Affine::IDENTITY, &grad1, None, &c,);
      } // ↓ in case step int conversion missed the last sliver
      let rad0_last = (skip_beg_deg + f64::from(steps_left) * precision_degps).to_radians();
      if rad0_last < r1end_rad { //println!("{rad0_last:.4}<{r1end_rad:.4} step_last {:.1}°<{:.1}° end",rad0_last.to_degrees(),r1end);
        let c = CircleSegment::new((cx,cy), r0,0.   ,  rad0_last,r1end_rad - rad0_last).outer_arc();
        let stroke_c = get_stroke_end(w_delta_avg);
        // scene.stroke(&stroke_c, Affine::IDENTITY, &grad1, None, &c,); // use col_avg? though grad should cover
        scene.stroke(&stroke_c, Affine::IDENTITY, &css::ORANGE, None, &c,); // for testing
      }


      let dash_off:f64 = 10.;
      // let dash_iter = [20.,20.,10.,10.]; //todo: reject negative numbers
      let dash_iter = [20.,20.]; //todo: reject negative numbers
      let deg_len = 2. * f64c::PI * r0 / 360.; //2π*100/360 = 1.74
      let rad_len = 2. * f64c::PI * r0 / 360.0_f64.to_radians(); //2π*100/6.28 = 100
      let dash_iter_len:f64 = dash_iter.iter().sum::<f64>(); //5
      let dash_iter_len_deg:f64 = dash_iter_len / deg_len; //° 2.87356
      let dash_off_deg = dash_off / deg_len;
      let dash_off_rad = dash_off_deg.to_radians();
      let dash_iter_len_rad = dash_iter_len_deg.to_radians(); //0.05015
      let dash_iter_rad = dash_iter.iter().map(|w| w / rad_len).collect::<Vec<f64>>();
      //3/(3+2) * 2.87 = 3/dash_iter_len * (dash_iter_len / deg_len) = 3 / deg_len
      //1.044

      // Draw pos-gradwidth segment separately without the extra iterator
      let c = CircleSegment::new((cx,cy), r0,0.   ,  r2beg_rad + rad_delta, skip_beg_rad).outer_arc();
      // let stroke_c = get_stroke_end(w2px);
      let stroke_c = get_stroke_end(w2px).with_dashes(dash_off,dash_iter);
      // scene.stroke(&stroke_c, Affine::IDENTITY, &col_end, None, &c,);
      scene.stroke(&stroke_c, Affine::IDENTITY, &css::WHEAT, None, &c,); // for testing

      // DEBUG copy shifted right
      // let cyy = cy;
      // let cxx = cx+5.;
      // let grad2cc_p0 = ( cxx + r0*f64::cos( r2beg_rad                      ) , cyy + r0*f64::sin( r2beg_rad                      ) );
      // let grad2cc_p1 = ( cxx + r0*f64::cos((r2beg + deg_delta).to_radians()) , cyy + r0*f64::sin((r2beg + deg_delta).to_radians()) );
      // let grad2cc = Gradient::new_linear(grad2cc_p0, grad2cc_p1).with_stops([col_avg    ,col_end]);
      // let c = CircleSegment::new((cxx,cyy), r0,r0,  r2beg_rad,arc_len_f.to_radians()).outer_arc();
      // let stroke_c = get_stroke_end(w2).with_dashes(dash_off,dash_iter);
      // scene.stroke(&stroke_c, Affine::IDENTITY, &grad2cc, None, &c,);

      // DEBUG copy smaller (including dashes, should perfectly align)
      let cyy = cy;
      let cxx = cx;
      let r00 = r0 - w2 - 1.;
      let grad2cc_p0 = ( cxx + r00*f64::cos( r2beg_rad                      ) , cyy + r00*f64::sin( r2beg_rad                      ) );
      let grad2cc_p1 = ( cxx + r00*f64::cos((r2beg + deg_delta).to_radians()) , cyy + r00*f64::sin((r2beg + deg_delta).to_radians()) );
      let grad2cc = Gradient::new_linear(grad2cc_p0, grad2cc_p1).with_stops([col_avg    ,col_end]);
      let c = CircleSegment::new((cxx,cyy), r00,r00,  r2beg_rad,arc_len_f.to_radians()).outer_arc();
      let stroke_c = get_stroke_end(w2).with_dashes(dash_off*r00/r0,dash_iter.iter().map(|w| w*r00/r0).collect::<Vec<f64>>());
      scene.stroke(&stroke_c, Affine::IDENTITY, &grad2cc, None, &c,);

      let sign2 = if w2 > wavg { 1.} else if w2 < wavg {-1.} else {0.}; //(from avg) ↑ if bigger, ↓ if smaller
      for i in 0..steps_left { let r = f64::from(i);
        let seg0 = r * precision_radps; // segment beginning in our arc coords (arc start = 0)
        let rad0 = r2beg_rad + seg0;
        let rad1 = rad0 + precision_radps; // todo debug only
        // let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,precision_radps+gap_correct).outer_arc(); //arc bugs with gaps
        // let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,precision_radps);
        //                          center  rout/in    ∠start ∠sweep
        let cw = wavg + sign2 * r * w_step_left;
        // let stroke_c = get_stroke_end(cw);
        let stroke_c = get_stroke_end(w2); // todo: same width for comparison with a reference

          // let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,precision_radps);
          // scene.stroke(&stroke_c, Affine::IDENTITY, &grad2, None, &c,);

        let seg_off = dash_off_rad % dash_iter_len_rad;
        let seg_beg = (dash_off_rad + seg0) % dash_iter_len_rad;// cut off arc that fit into the previous dash set, so this is our segment beginning in the coordinate system of a dash set (set's begin = 0)
        let seg_end = seg_beg + precision_radps;
        let mut is_drawn = true;
        let mut d_beg = 0.; // length up to the beginning of this dash = ∑ of all previous dash lens

        //     ──────  —————  outer line dash pattern (todo: what if our line is bigger than 1 pattern? like here)
        //        ┌─────┐     our   line (always starts later due to "rad0 % dash_iter_len_rad")
        // seg_beg┘     └seg_end
        //         ↑↑  ↑ draw, overlaps with   active
        //           ↑↑  skip, overlaps with inactive
        if i == 0 {println!("\n\n—————————————————————————————————————————————————————————————————————————————————")};
        for dash_i in &dash_iter_rad {
          if is_drawn { // ignore inactive dashes
            let d_end = d_beg + dash_i;
            let draw_beg = d_beg.max(seg_beg).min(d_end); // start at dash begin, → to segment begin, but not past dash end
            let draw_end = d_end.min(seg_end).max(d_beg); // start at dash end  , ← to segment end  , but not past dash beg
            let draw_len = draw_end - draw_beg;
            // if rad0      <=       d_end
              // &&    seg_end >= d_beg  { // our segment overlaps with this dash
            println!(
              "{}abs {: >4.1}° → {: >4.1}° Δ{: >3.1}° off {: >3.1}°¦{: >3.1}°\
              │ rel {: >4.1}° → {: >4.1}° Δ{: >3.1}°\
              │ dash {: >4.1}° → {: >4.1}° Δ{: >4.1}°\
              │ draw {: >4.1}° → {: >4.1}° ⇒ {: >3.1}° "
              ,if draw_len>0.{"✓ "}else{"  "}
              ,rad0    .to_degrees(),rad1    .to_degrees(),(rad1-rad0).to_degrees(), dash_off_deg, seg_off
              ,seg_beg .to_degrees(),seg_end .to_degrees(),(seg_end - seg_beg).to_degrees()
              ,d_beg   .to_degrees(),d_end   .to_degrees(),dash_i.to_degrees()
              ,draw_beg.to_degrees(),draw_end.to_degrees(),draw_len.to_degrees()
              );
            if draw_len > 0. {
              let c = CircleSegment::new((cx,cy), r0,r0   ,rad0,draw_len);
              scene.stroke(&stroke_c, Affine::IDENTITY, &grad2, None, &c,);
            }
          } //else {println!("   inactive ({: >4.1}°)",dash_i.to_degrees());}
          d_beg += dash_i;
          is_drawn = !is_drawn;
        }
      } // ↓ in case step int conversion missed the last sliver
      let rad0_last = (r2beg + f64::from(steps_left) * precision_degps).to_radians();
      if rad0_last < skip_beg_rad { //println!("{rad0_last:.4} < {skip_beg_rad:.4} step_last {:.1}°<{:.1}° end",rad0_last.to_degrees(),skip_beg_deg);
        let c = CircleSegment::new((cx,cy), r0,0.   ,  rad0_last,skip_beg_rad - rad0_last).outer_arc();
        let stroke_c = get_stroke_end(w2px);
        // scene.stroke(&stroke_c, Affine::IDENTITY, &grad2, None, &c,); // use col_avg? though grad should cover
        scene.stroke(&stroke_c, Affine::IDENTITY, &css::WHEAT, None, &c,); // for testing
      }

      // Draw debug circles showing where each gradient begins/ends
      let pstr = get_stroke_end(1.); // starting point bigger than the ending, angle to differentiate two curves
      let g1beg = Ellipse::new(grad1_p0, ( 2.,15.+5.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g1beg,);
      let g1end = Ellipse::new(grad1_p1, ( 2.,15.+0.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g1end,);
      let g2beg = Ellipse::new(grad2_p0, ( 2.,15.+5.),  99.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g2beg,);
      let g2end = Ellipse::new(grad2_p1, ( 2.,15.+0.),  99.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g2end,);

    }
  }
}
