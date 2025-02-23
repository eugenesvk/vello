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
  use crate::SceneParams;
  use kurbo::RoundedRect;
  use vello::kurbo::{Affine, BezPath, Cap, Arc,Circle,CircleSegment,Ellipse,Line,Join, PathEl, Point, Rect, Shape, Stroke, Vec2,};
  use vello::peniko::color::{palette, AlphaColor, Lch, palette::css};
  use vello::peniko::*;
  use vello::*;
  use std::f64::consts as f64c;
  const epsi:f64 = 0.00000000001; // to counter some float precision errors

  fn get_stroke    (width:f64) -> Stroke {Stroke::new(width).with_start_cap(Cap::Butt).with_end_cap(Cap::Round).with_join(Join::Bevel)}
  fn get_stroke_end(width:f64) -> Stroke {Stroke::new(width).with_start_cap(Cap::Butt).with_end_cap(Cap::Butt ).with_join(Join::Bevel)}
  fn get_stroke1   (width:f64) -> Stroke {Stroke::new(width).with_start_cap(Cap::Butt).with_end_cap(Cap::Square).with_join(Join::Round)}

  // Two curves joined by converging colors and widths, with an example of bottom/top Arcs:
    // Curve is split in 2 parts: 1st near joins (≝~join) (end for the bottom arc, beg for the top) 2nd the rest
    // ~Joins part is styled manually, the rest is styled as a regular line
    // Color transition towards average:
      // Regular "gradient" works, start at ~join outer edge with line's base color, end at the joint with an average color between lines
      // TODO: if main lines have gradients themselves, this doesn't work, and not sure what the proper logic would be to glue arbitrary gradient styles together (maybe add handling of a few common cases?)
    // Widths transition towards average:
      // split each ~join into X steps based on user defined precision
      // for each step, get the width delta and draw a small arc with the base width + steps * w_delta
        // use Circle segment with both r1=r2 to avoid occlusion artifacts or make each arc end a bit futher to overalp
      // last step is different as it doesn't cover the same length (since not all lines can be split into a whole number of equal steps)
    // Dashes don't transition, not sure whether any logic that would introduce a 3rd var pattern is better:
      // (this just describes the approach of adding dashes to our step-by-step approach)
      // calculate dashed set [1,2,1,4] length in radians (1+2+1+4=8px == 1 rad at radius X)
      // for each step, get where each step begin lands in set coordinates by calculating division remainder of step begin in curve coordinates (+offset) by set length
      // within each set, iterate by dash item, and if it's an active/drawn item, determine how much of our step length fits there, and draw it
    // Dashes harder/failed approaches: combine vertical line with 1/2 circle; add dashed strokes; split by path and try to manually check which should be drawn;  change stroke width in the last segments if they are covered by the 1/2 circle (not possible???? how to calculate a match? check if we can draw a semicircle and check if overlap); draw previous line;  draw splits with a different width and gradient

  // TODO
    // convert circle segments to Arcs directly and overlap (except for the last segment) to avoid conflaction artifacts
    // test if step length > dash set length (with very low precision)
    // reject negative numbers on accepted dash iterator
    // ?? convert rads to degree floats to avoid small errors on adding dashes?
    // + calculate the remainder from iterative approach and use it as a (-) offset to the main curve
  // ??? update offset algo to find index to the dash that matches offset ???
    //  let mut dash_ix = 0;
    //  let mut dash_remaining = dashes[dash_ix] - dash_offset;
    //  let mut is_active = true;
    //  // Find place in dashes array for initial offset.
    //  while dash_remaining < 0.0 {
    //      dash_ix = (dash_ix + 1) % dashes.len();
    //      dash_remaining += dashes[dash_ix];
    //      is_active = !is_active;
    //  }

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
      // Size
      let arc_len_deg:f64 = 180.;
      let precision_deg_per_step:f64 = 0.5; //rad 0.00873
      let steps_f = arc_len_deg / precision_deg_per_step; //360
      let steps_i   = steps_f as i32;
      // Gradient / size convergence bounds
      let f_delta = 0.337; // start changing width for the first/last % only
      let delta_deg    	= arc_len_deg *       f_delta ; let rad_delta    = delta_deg   .to_radians();
      let skip_beg_deg 	= arc_len_deg * (1. - f_delta); let skip_beg_rad = skip_beg_deg.to_radians();
      let steps_delta_f	= steps_f   * f_delta; //36
      let skip_end     	= steps_delta_f as i32; //36
      let skip_beg     	= steps_i - skip_end; //324

      // ((w1 + sign1 * w_step * r) * dpi).round() / dpi; //transition shouldn't be pixel-stepped!
      let gap:f64 = 0.; // doesn't seem to 0.0001 affect anything with corrected ending style to Bevel
      let r1beg:f64 = 0.             	; let r1beg_rad = r1beg.to_radians(); //→
      let r1end = r1beg + arc_len_deg	; let r1end_rad = r1end.to_radians();
      let r2beg = r1end + gap        	; let r2beg_rad = r2beg.to_radians();
      let r2end = r2beg + arc_len_deg	; let r2end_rad = r2end.to_radians();
      let col_beg = css::LIME;
      let col_end = css::RED;
      let col_avg = col_beg.lerp(col_end,0.5,Default::default());

      // Segment 1: ~join part is 2nd (at the end)
      // Draw pre-gradwidth segment separately without the extra iterator
      let c = CircleSegment::new((cx,cy), r0,0.   ,  r1beg_rad,skip_beg_rad).outer_arc();
      let stroke_c = get_stroke_end(w1px);
      // scene.stroke(&stroke_c, Affine::IDENTITY, &col_beg, None, &c,);
      scene.stroke(&stroke_c, Affine::IDENTITY, &css::ORANGE, None, &c,); // for testing

      let steps_delta_f     	= delta_deg / precision_deg_per_step; // 180°*33% / 0.5 = 118.8
      let steps_delta_i     	= steps_delta_f as i32; let steps_delta_if = f64::from(steps_delta_i); // 118 whole steps for later iteration and drawing by small step
      let delta_covered_deg 	=                  steps_delta_if  * precision_deg_per_step;
      // let steps_delta_rem	=  steps_delta_f - steps_delta_if; //0.8 steps not covered by whole
      let delta_rem_deg     	= delta_deg - delta_covered_deg; let delta_rem_rad = delta_rem_deg.to_radians();

      let w_per_step_i:f64	= w_delta_avg / f64::from(steps_delta_i); //to reach average

      let sign1 = if w1 > wavg {-1.} else if w1 < wavg {1.} else {0.}; //(to avg) ↓ if bigger, ↑ if smaller
      for i in 0..steps_delta_i { let r = f64::from(i);
        let rad0 = r2beg_rad + r * precision_rad_per_step;
        let c = CircleSegment::new((cx,cy), r0,0.   ,  rad0,precision_rad_per_step).outer_arc();
        //                          center  rout/in    ∠start ∠sweep
        // if i == 0 {println!("\n\n———————————————————————————————")};
        // println!("i={i}/{steps_delta_i}  r={r:.1} r1beg={r1beg:.1} r*prec={:.1} deg={skip_beg_deg:.1} beg={:.1} end={:.1}",r * precision_deg_per_step
        //   ,       skip_beg_deg + r * precision_deg_per_step, skip_beg_deg + r * precision_deg_per_step + precision_deg_per_step);
        let cw = w1 + sign1 * r * w_per_step_i;
        let stroke_c = if i >=  (steps_delta_i - 3) {get_stroke_end(cw) //last steps no round ends
        } else {get_stroke(cw)}; // round ends to fix the outer arc artifacts when joining in rectangles without curves
        scene.stroke(&stroke_c, Affine::IDENTITY, &grad1, None, &c,);
      } // ↓ in case step int conversion missed the last sliver
      let rad0_last = (skip_beg_deg + f64::from(steps_delta_i) * precision_deg_per_step).to_radians();
      if rad0_last < r1end_rad { //println!("{rad0_last:.4}<{r1end_rad:.4} step_last {:.1}°<{:.1}° end",rad0_last.to_degrees(),r1end);
        let c = CircleSegment::new((cx,cy), r0,0.   ,  rad0_last,r1end_rad - rad0_last).outer_arc();
        let stroke_c = get_stroke_end(w_delta_avg);
        // scene.stroke(&stroke_c, Affine::IDENTITY, &grad1, None, &c,); // use col_avg? though grad should cover
        scene.stroke(&stroke_c, Affine::IDENTITY, &css::ORANGE, None, &c,); // for testing
      }

      // Position
      let cx = 900.; let cy = 200.; let r0 = 95.5; //600 circum len 300 half
      let deg_len = 2. * f64c::PI * r0 / 360.; //2π*100/360 = 1.74

      let dash_off_deg = 30.1; let dash_iter_deg = [30.1,40.];
      let dash_iter = dash_iter_deg.iter().map(|w|w*deg_len).collect::<Vec<f64>>();
      let dbg = 0;

      let wpx = w2px;
      ddd(scene, (cx,cy),r0, r2beg_rad, JoinWhere::Beg,
        col_avg,col_end,
        wpx,wavg,
        steps_delta_i,w_per_step_i,precision_rad_per_step,
        delta_rem_rad,
        dash_off_deg,dash_iter_deg.to_vec(),
        rad_delta,
        skip_beg_rad,
        delta_deg,arc_len_deg, dbg,
      );
      ddd_debug(scene, (cx,cy),r0, r2beg_rad, arc_len_deg, col_avg,col_end, wpx, delta_deg, dash_off_deg,dash_iter_deg.to_vec());

      // Draw debug circles showing where each gradient begins/ends
      let pstr = get_stroke_end(1.); // starting point bigger than the ending, angle to differentiate two curves
      let g1beg = Ellipse::new(grad1_p0, ( 2.,15.+5.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g1beg,);
      let g1end = Ellipse::new(grad1_p1, ( 2.,15.+0.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g1end,);
      let g2beg = Ellipse::new(grad2_p0, ( 2.,15.+5.),  99.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g2beg,);
      let g2end = Ellipse::new(grad2_p1, ( 2.,15.+0.),  99.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g2end,);

    }
  }
  pub enum JoinWhere{Beg,End,}
  use std::fmt::Display;
  impl Display for JoinWhere {fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { match self {
    JoinWhere::Beg	=> write!(f, "╍╍——"),
    JoinWhere::End	=> write!(f, "——╍╍")}   }
  }

  // DEBUG copy smaller (including dashes, should perfectly align as dash length/offsets are adjusted per difference in size)
  pub fn ddd_debug(scene:&mut Scene, center:impl Into<Point>,r0:f64, arc_beg:f64, arc_len_deg:f64,
    col_avg:Color,col_end:Color,
    wpx:f64,
    delta_deg:f64,
    dash_off_deg:f64,dash_iter_deg:Vec<f64>,
    ) {
    let c = center.into(); let (cxx,cyy) = (c.x,c.y); let r00 = r0 - wpx;
    let deg_len  	= 2. * f64c::PI * r0 / 360.0                 ; //2π*100/360 = 1.74
    let rad_len  	= 2. * f64c::PI * r0 / 360.0_f64.to_radians(); //2π*100/6.28 = 100
    let dash_off 	= dash_off_deg * deg_len;
    let dash_iter	= dash_iter_deg.iter().map(|w|w*deg_len).collect::<Vec<f64>>();

    let grad2cc_p0 = ( cxx + r00*f64::cos(arc_beg                         ) , cyy + r00*f64::sin(arc_beg                         ) );
    let grad2cc_p1 = ( cxx + r00*f64::cos(arc_beg + delta_deg.to_radians()) , cyy + r00*f64::sin(arc_beg + delta_deg.to_radians()) );
    let grad2cc = Gradient::new_linear(grad2cc_p0, grad2cc_p1).with_stops([col_avg    ,col_end]);
    let c = Arc::new((cxx,cyy), (r00,r00),  arc_beg,arc_len_deg.to_radians(), 0.);
    let stroke_c = get_stroke_end(wpx).with_dashes(dash_off*r00/r0,dash_iter.iter().map(|w| w*r00/r0).collect::<Vec<f64>>());
    scene.stroke(&stroke_c, Affine::IDENTITY, &grad2cc, None, &c,);

  }
  pub fn ddd(scene:&mut Scene, center:impl Into<Point>,r0:f64,  arc_beg:f64, jn:JoinWhere,
    col_avg:Color,col_end:Color,
    wpx:f64,wavg:f64,
    steps_delta_i:i32,w_per_step_i:f64,precision_rad_per_step:f64,
    delta_rem_rad:f64,
    dash_off_deg:f64,dash_iter_deg:Vec<f64>,
    rad_delta:f64,
    skip_beg_rad:f64,
    delta_deg:f64,
    arc_len_deg:f64,     dbg:u8,
    ) {
      let mut is_vis_draw  = false; // whether a visible dash is drawing to clamp its first partial draw to the next. Switches to off when an invisible dash is "drawing"
      let mut dash_partial = 0.; // use as dash offset for the next segment to hide the partially drawn part
      let mut carry_over:f64 = 0.; // if Δstep covers 2 dash segments, the 1st one will store the remainer it didn't cover here for the 2nd to pick it up

      let sign = match jn {
        JoinWhere::Beg	=> if wpx > wavg { 1.} else if wpx < wavg {-1.} else {0.}, //(from avg) ↑ if bigger, ↓ if smaller
        JoinWhere::End	=> if wpx > wavg {-1.} else if wpx < wavg { 1.} else {0.}, //(to   avg) ↓ if bigger, ↑ if smaller
      };
      let c = center.into(); let (cx,cy) = (c.x,c.y);

      let deg_len = 2. * f64c::PI * r0 / 360.0                 ; //2π*100/360 = 1.74
      let rad_len = 2. * f64c::PI * r0 / 360.0_f64.to_radians(); //2π*100/6.28 = 100
      let dash_off         	= dash_off_deg * deg_len;
      let dash_iter        	= dash_iter_deg.iter().map(|w|w*deg_len).collect::<Vec<f64>>();
      let dash_iter_len_px 	= dash_iter.iter().sum::<f64>(); //5
      let dash_iter_len_deg	= dash_iter_len_px / deg_len; let dash_iter_len_rad = dash_iter_len_deg.to_radians(); //2.87356° 0.05015
      let dash_off_deg     	= dash_off / deg_len; let dash_off_rad = dash_off_deg.to_radians();
      let dash_iter_rad = dash_iter.iter().map(|w| w / rad_len).collect::<Vec<f64>>();

      let (grad_p0,grad_p1) = match jn {
        JoinWhere::Beg	=> {
          let grad_beg_p0 = ( cx + r0*f64::cos(arc_beg                           ) , cy + r0*f64::sin(arc_beg                         ) );
          let grad_beg_p1 = ( cx + r0*f64::cos(arc_beg + delta_deg.to_radians()  ) , cy + r0*f64::sin(arc_beg + delta_deg.to_radians()) );
          (grad_beg_p0,grad_beg_p1)},
        JoinWhere::End	=> {
          let grad_end_p0 = ( cx + r0*f64::cos(skip_beg_rad                      ) , cy + r0*f64::sin(skip_beg_rad));
          let grad_end_p1 = ( cx + r0*f64::cos(arc_beg + arc_len_deg.to_radians()) , cy + r0*f64::sin(arc_beg + arc_len_deg.to_radians()));
          (grad_end_p0,grad_end_p1)},
      };
      let gap:f64 = 0.; // doesn't seem to 0.0001 affect anything with corrected ending style to Bevel
      let r1beg:f64 = 0.             	; let r1beg_rad = r1beg.to_radians(); //→
      let r1end = r1beg + arc_len_deg	; let r1end_rad = r1end.to_radians();
      let r2beg = r1end + gap        	; let r2beg_rad = r2beg.to_radians();

      let grad = Gradient::new_linear(grad_p0,grad_p1).with_stops([col_avg,col_end]);



      if let JoinWhere::End = jn {// Segment 1: ~join part is 2nd (at the end)
        // Draw pre-gradwidth segment separately without the extra iterator
        let c = Arc::new((cx,cy), (r0,r0)   ,  arc_beg,skip_beg_rad, 0.);
        let stroke_c = get_stroke_end(wpx);
        // scene.stroke(&stroke_c, Affine::IDENTITY, &col_beg, None, &c,);
        scene.stroke(&stroke_c, Affine::IDENTITY, &css::ORANGE, None, &c,); // for testing
      }

      // TODO: force precision to be so that one step is never bigger than a dash set length, otherwise would need to repeat the full dash parsing logic here
        // or maybe have a better loop logic that can handle it?
      let is_extra_step = delta_rem_rad > 0.;
      let steps_delta_xt = if is_extra_step {steps_delta_i} else {steps_delta_i-1};
      for i in 0..=steps_delta_xt { let r = f64::from(i); let is_last = i == steps_delta_xt;
        // NB! last step needs special handling since it's a fractional one, so not full "precision length"!
        let step_width = if is_last && is_extra_step	{delta_rem_rad
        } else                                      	{precision_rad_per_step};
        let seg0 = (r - 1.) * precision_rad_per_step + step_width; // segment beg in our arc coords (arc start = 0), replace last regular width with a custom step_width
        let rad0 = arc_beg + seg0;
        let rad1 = rad0 + step_width; // todo debug only
        // let c = Arc::new((cx,cy), (r0,r0) ,  rad0,step_width+gap_correct, 0.); //arc bugs with gaps
        // let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,step_width); // alt fix
        //                          center  rout/in    ∠start ∠sweep
        let cw = if is_last	{wpx
        } else             	{wavg + sign * r * w_per_step_i};
        // let stroke_c = get_stroke_end(cw);
        let stroke_c = get_stroke_end(wpx); // todo: same width for comparison with a reference

          // let c = Arc::new((cx,cy), (r0,r0)   ,  rad0,step_width, 0.);
          // scene.stroke(&stroke_c, Affine::IDENTITY, &grad, None, &c,);

        let seg_off =  dash_off_rad         % dash_iter_len_rad;
        let seg_be_ = (dash_off_rad + seg0) % dash_iter_len_rad;// cut off arc that fit into the previous dash set, so this is our segment beginning in the coordinate system of a dash set (set's begin = 0)
        let (seg_beg,seg_count) = if (dash_iter_len_rad - seg_be_).abs() < epsi { // segment starts at dash set end's, so shift it to next dash set's begin (float imprecision prevents clean division)
          (0.     , (dash_off_rad + seg0).div_euclid(dash_iter_len_rad) + 1.+1.)
        } else {
          (seg_be_, (dash_off_rad + seg0).div_euclid(dash_iter_len_rad) + 1.   )
        };
        let seg_count = (dash_off_rad + seg0).div_euclid(dash_iter_len_rad) + 1.;
        let seg_end = seg_beg + step_width;
        let mut d_beg = 0.; // length up to the beginning of this dash = ∑ of all previous dash lens

        //     ──────  —————  outer line dash pattern (todo: what if our line is bigger than 1 pattern? like here)
        //        ┌─────┐     our   line (always starts later due to "rad0 % dash_iter_len_rad")
        // seg_beg┘     └seg_end
        //         ↑↑  ↑ draw, overlaps with   active
        //           ↑↑  skip, overlaps with inactive
        // if i == 0 {println!("\n\n—————Σⁱ={steps_delta_xt: >3}——╍╍№{} Σ{dash_iter_len_deg: >4.1}° off{dash_off_deg: >4.1}° {dash_iter_deg:?}°╍╍——beg {r2beg: >4.1}° Δ{delta_covered_deg: >4.1}° + {delta_rem_deg: >4.1}° rem = Δ{delta_deg: >4.1}°————————————————————————————————————————————"
          // ,dash_iter.len())};
        let mut dr = 0; // track dash 🗘
        let mut step_covered = step_width; // track Σ dash_iter_len_rad covering each Δstep
        while step_covered > 0.  {
          dr += 1;
          step_covered -= dash_iter_len_rad;
        let mut j = 0;
        let mut is_visible = false;
        let mut prev_draw_len:f64 = 0.; // store a sum of previously drawn dashes (vis+invis) so that if this dash doesn't cover the full Δstep or Σdash_len, we can see which part of it was covered before and which should go as Δover to the next step
        // if seg_count == 1. {
        for dash_i in &dash_iter_rad {
          is_visible = !is_visible;
          j += 1;
          // if is_visible && j == 3 {
          // if seg_end < d_beg {
          //   println!("{: >4.1}° {: >4.1}° → {: >4.1}°",seg_end.to_degrees(),d_beg.to_degrees(),(d_beg + dash_i).to_degrees());
          //   break;
          // } // our segment has been fully covered, no need to continue
          let mut dbgprint = false;
          if is_visible {
            let d_end = d_beg + dash_i;
            let draw_beg = d_beg.max(seg_beg).min(d_end); // start at dash begin, → to segment begin, but not past dash end
            let draw_end = d_end.min(seg_end).max(d_beg); // start at dash end  , ← to segment end  , but not past dash beg
            let draw_len = draw_end - draw_beg;

            if carry_over > 0. { // draw leftovers from the previous dash
              prev_draw_len += carry_over;
              // if is_vis_draw {println!("!!! leftovers from a previous dash should always 1st, but something else drew")}; //todo warn
              is_vis_draw = true; // start drawing @ the end of prev ↓ step
              let c = Arc::new((cx,cy), (r0,r0)   ,rad0 - carry_over,carry_over, 0.);
              // scene.stroke(&stroke_c, Affine::IDENTITY, &grad    , None, &c,);
              scene.stroke(&stroke_c, Affine::IDENTITY, css::MAGENTA , None, &c,); // todo: replace test with ↑
              carry_over = 0.;
              dbgprint = true; // todo: remove
              // println!("{i} drawing Δover {: >2.1} @ {: >3.2} = ({: >2.1}-{: >2.1}-Δ{: >2.1}) → {: >3.2}",carry_over.to_degrees()
                // ,(rad1 - step_width - carry_over).to_degrees(),rad1.to_degrees(),step_width.to_degrees(),carry_over.to_degrees()
                // ,(rad1 - step_width).to_degrees());
            }
            if draw_len > 0.0 { // 1st draw starts @ seg end to attach to the next draw in case of partials
              prev_draw_len += draw_len;
              let c = if is_vis_draw   {Arc::new((cx,cy), (r0,r0)   ,rad0         ,draw_len, 0.)
              } else {is_vis_draw=true; Arc::new((cx,cy), (r0,r0)   ,rad1-draw_len,draw_len, 0.)};
              if is_last	{scene.stroke(&stroke_c, Affine::IDENTITY, &css::LIME, None, &c,);
              } else    	{scene.stroke(&stroke_c, Affine::IDENTITY, &grad    , None, &c,);}
              // scene.stroke(&stroke_c, Affine::IDENTITY, &grad, None, &c,); // todo: replace ↑ lime test with ←
              if is_last && draw_len < *dash_i - epsi { // drawn something, but not the full visible dash
                let part_len = draw_end - d_beg; //how much of an existing dash is covered by all draws, incl. last
                dash_partial = (d_beg + part_len) * rad_len; // add all prior dash segments within a set
                // println!("! last +visible +draw ╍w {: >.4}° − {: >.4}° actual = {: >.4}° {: >.4}°part_len partial ({:.1}px) rad1 {:.3}°"
                //   ,dash_i.to_degrees(),draw_len.to_degrees(), (d_end.min(seg_end) - d_beg).to_degrees(),part_len.to_degrees(), dash_partial,rad1.to_degrees());
              }
            } else {is_vis_draw=false;}
            // if rad0       <=       d_end
            // &&    seg_end >= d_beg  { // (alt check) our segment overlaps with this dash
            if dbg>0 && (dbgprint || i == 0 || is_last || (78<= i && i <=83)) {println!( //👁👀👓  seg={dash_off_deg: >3.1} % {dash_iter_len_deg: >3.1}
              "{}👀{}{i: >3} {} {j: >2}\
              │ +{: >4.1}={: >4.1}° ↷ {: >4.1}° Δ{: >3.1}° off {: >3.1}° \
              │№{seg_count: >2} {: >4.1}° ↷ {: >4.1}°\
              │ ╍ {: >4.1}° → {: >4.1}° Δ{: >4.1}°\
              │ 🖉 {: >4.1}° → {: >4.1}° ⇒ {: >3.1}° "
              ,if draw_len>0.{"🖉"}else{" "}, if is_last {"🛑"}else{" "}, if dr > 1 {format!("🗘{dr: >1}")}else{"  ".to_owned()}
              ,seg0.to_degrees()    ,rad0    .to_degrees(),rad1    .to_degrees(),(rad1    -    rad0).to_degrees(),seg_off.to_degrees()
              ,                      seg_beg .to_degrees(),seg_end .to_degrees()
              ,d_beg   .to_degrees(),d_end   .to_degrees(),dash_i.to_degrees()
              ,draw_beg.to_degrees(),draw_end.to_degrees(),draw_len.to_degrees()
              );}
          } else { //println!("   inactive ({: >4.1}°)",dash_i.to_degrees());
            let d_end = d_beg + dash_i;
            let draw_beg = d_beg.max(seg_beg).min(d_end); // start at dash begin, → to segment begin, but not past dash end
            let draw_end = d_end.min(seg_end).max(d_beg); // start at dash end  , ← to segment end  , but not past dash beg
            let draw_len = draw_end - draw_beg;
            if draw_len > 0.0 {is_vis_draw = false; prev_draw_len += draw_len;
              let space_available = step_width.min(dash_iter_len_rad) - prev_draw_len;
              if space_available > 0.00000000001 { // this+prev dashes didn't cover the full Δstep¦dash segment width (whichever is smaller, if dash segment fits in Δstep, then ), so should be drawn by the next visible dash
                carry_over = space_available;
                if is_last { // no next step, draw in this one
                  is_vis_draw = true;
                  let carry_over_r0 = rad0 + prev_draw_len;
                  let carry_over_r1 = (carry_over_r0 + carry_over).min(rad1) ;// up to our arc's end, the rest will be picked up by the next arc
                  let carry_over_delta = carry_over_r1 - carry_over_r0;
                  let c = Arc::new((cx,cy), (r0,r0)   ,carry_over_r0,carry_over_delta, 0.);
                  // scene.stroke(&stroke_c, Affine::IDENTITY, &grad    , None, &c,);
                  scene.stroke(&stroke_c, Affine::IDENTITY, css::MAGENTA , None, &c,); // todo: replace test with ↑
                  dash_partial = carry_over_delta * rad_len; carry_over = 0.;
                  // println!("last step - drawn next dash since it won't be handled later!");
                // } else {println!("  Δover {: >4.1}° = step_w {: >4.1}° - {: >4.1}° drawn",carry_over.to_degrees(),step_width.to_degrees(),draw_len.to_degrees());
                }
              }
            }
            if is_last { // last invisible dash also needs to signal its width to update offset of the next arc
              let part_len = draw_end - d_beg; //how much of an existing dash is covered by all draws, incl. last
              if   draw_len > 0. //drawn something… ↙some float rounding error
                && part_len < *dash_i - 0.00000000001 { //…but not the full invisible dash
                dash_partial = (d_beg + part_len) * rad_len; //≝draw_end add all prior dash segments within a set
                // println!("{}№{seg_count} last -visible +draw ╍beg {: >3.3}° draw_end {: >3.3}°≝partial ({:.1}px) Δ{: >3.3}° Δstep {: >3.3}° drawn │ ╍w {: >.2}° left {: >.2}° rad1 {:.3}°"
                //   ,if dash_partial > 0. {"✓"}else{"✗"}
                //   ,d_beg.to_degrees(),draw_end.to_degrees(),dash_partial,part_len.to_degrees(),draw_len.to_degrees()
                //   ,dash_i.to_degrees(),(dash_i-part_len).to_degrees(),rad1.to_degrees());
              }
            }
            if dbg>1 && (dbgprint || i == 0 || is_last || (78<= i && i <=83)) {println!( //👁👀👓  seg={dash_off_deg: >3.1} % {dash_iter_len_deg: >3.1}
              "{}👓{}{i: >3} {} {j: >2}\
              │ +{: >4.1}={: >4.1}° ↷ {: >4.1}° Δ{: >3.1}° off {: >3.1}° \
              │№{seg_count: >2} {: >4.1}° ↷ {: >4.1}°\
              │ ╍ {: >4.1}° → {: >4.1}° Δ{: >4.1}°\
              │ 🖉 {: >4.1}° → {: >4.1}° ⇒ {: >3.1}° "
             ,if draw_len>0.{"🖉"}else{" "}, if is_last {"🛑"}else{" "}, if dr > 1 {format!("🗘{dr: >1}")}else{"  ".to_owned()}
              ,seg0.to_degrees()    ,rad0    .to_degrees(),rad1    .to_degrees(),(rad1    -    rad0).to_degrees(),seg_off.to_degrees()
              ,                      seg_beg .to_degrees(),seg_end .to_degrees()
              ,d_beg   .to_degrees(),d_end   .to_degrees(),dash_i.to_degrees()
              ,draw_beg.to_degrees(),draw_end.to_degrees(),draw_len.to_degrees()
              );}
          }
          d_beg += dash_i;
        }
        }
      }

      if let JoinWhere::Beg = jn { // Segment 2: ~join part is 1st (at the start)
        // Draws pos-gradwidth segment separately without the extra iterator, including leftovers from whole steps not covering the full range
        // println!("2nd curve@end: {: >3.0}° + Δ{: >3.0}° ⇒ {: >3.0}° + {: >3.0} skip_beg ╍part={: >3.0}° ({dash_partial: >3.0}px)"
        //   ,arc_beg.to_degrees(),rad_delta.to_degrees(),(arc_beg + rad_delta).to_degrees(),skip_beg_rad.to_degrees(),(dash_partial/rad_len).to_degrees());
        let c = Arc::new((cx,cy), (r0,r0) , arc_beg + rad_delta, skip_beg_rad, 0.);
        // let stroke_c = get_stroke_end(wpx); // todo↓ make dashes optional
        let stroke_c = get_stroke_end(wpx).with_dashes(dash_partial,dash_iter); // use remainder from the previous segment so that the total matches the style as though it were drawn in one step
        // scene.stroke(&stroke_c, Affine::IDENTITY, &col_end, None, &c,);
        scene.stroke(&stroke_c, Affine::IDENTITY, &css::WHEAT, None, &c,); // for testing
      }
    }
}