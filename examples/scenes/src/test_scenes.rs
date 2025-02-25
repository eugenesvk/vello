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
  use crate::color::Srgb;
  use std::borrow::Borrow;
  use std::f64::consts as f64c;

  const epsi:f64 = 0.00000000001; // to counter some float precision errors
  // buggy imprecise, doesn't handle big numbers, but doesn't matter much here
  pub fn f_round(f:f64, digits:i32) -> f64 {let dig_pow=10_f64.powi(digits); (f * dig_pow).round() / dig_pow}

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

  // TODO:
    // don't extend last visible dash not to bleed into the invisible one
    // test if step length > dash set length (with very low precision)
    // maybe add a min average gap between two lines so that if first line ends with a partial inactive gap, the 2nd doesn't immediateely start with a visible dash, but + offset (unless it's too big for the 2d line, thus min average? or just min)
    // + convert circle segments to Arcs directly
      // add overlap (except for the last segment) to avoid conflaction artifacts
    // test if step length > dash set length (with very low precision)
    // reject negative numbers on accepted dash iterator
    // ?? convert rads to degree floats to avoid small errors on adding dashes?
    // + calculate the remainder from iterative approach and use it as a (-) offset to the main curve
    // + make gradient sweeps for better precision with arcs instead of linear
      // + for non-stepped lines, can use multiple colors with steps as well
  // Kurbo precision bug leading to artifacts comparing to a reference dashed circle (which is incorrect, ↓ block 4 from the bottom has incorrect shape vs other blocks)
    // let dash_off_deg = 0.; let dash_iter_deg = [10.,10.]; ← translates to >10e-6 numbers, while kurbo precision limit is e-6
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

  pub enum JoinWhere{Beg,End,}
  use std::fmt::Display;
  impl Display for JoinWhere {fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { match self {
    JoinWhere::Beg	=> write!(f, "╍╍——"),
    JoinWhere::End	=> write!(f, "——╍╍")}   }
  }

  pub(super) fn stroke_styles(transform: Affine) -> impl FnMut(&mut Scene, &mut SceneParams<'_>) {
    use PathEl::*;
    move |scene, params| {
      let dpi = 1.5;
      // Size
      let arc_len_deg:f64 = 180.; let precision_deg_per_step = 2.5; //rad 0.00873
      let delta_transit = 1.0; // start changing width for the first/last % only
      // Line width
      let w1:f64 = 20.; let w2:f64 =  4.;
      // Position
      let cx = 0.; let cy = 400.; let r0 = 395.5; //600 circum len 300 half
      let r1beg:f64 = 0.             	; let r1beg_rad = r1beg.to_radians(); //→
      let r2beg = r1beg + arc_len_deg	; let r2beg_rad = r2beg.to_radians();
      // Color
      let col_beg = css::LIME;let col_end = css::RED;
      // Dashes
      let mut dashes: Vec<(f64,Option<[f64;2]>, u8)> = vec![];
      // dash     offset, iter            ,dbg
      dashes.push( ( 0. , None            , 0) );
      dashes.push( ( 0. , None            , 1) );
      dashes.push( (30.1, Some([30.1,40.]), 1) );
      dashes.push( (15.1, Some([10.1,10.]), 1) );
      dashes.push( (12.1, Some([20.1,20.]), 1) );
      dashes.push( (11.1, Some([30.1,10.]), 1) );
      dashes.push( (0.  , Some([30.1,10.]), 1) );

      let end = dashes.len();
      for i in 0..end { let f = f64::from(i as u32);
        let cx = 20. + (1. + (f     %      5.) * 2. )*(r0 + w1.max(w2));
        let cy = 20. + (1. +  f.div_euclid(5.) * 2. )*(r0 + w1.max(w2)); // 5 circles in a row
        ddd(scene, (cx,cy),r0, r1beg_rad, arc_len_deg,  JoinWhere::End,delta_transit,
          col_beg,col_end,  w1,w2, dpi,  precision_deg_per_step,  dashes[i].0,dashes[i].1, dashes[i].2,);
        ddd(scene, (cx,cy),r0, r2beg_rad, arc_len_deg,  JoinWhere::Beg,delta_transit,
          col_beg,col_end,  w1,w2, dpi,  precision_deg_per_step,  dashes[i].0,dashes[i].1, dashes[i].2,);
      }

      // Draw debug circles showing where each gradient begins/ends
      // let pstr = get_stroke_end(1.); // starting point bigger than the ending, angle to differentiate two curves
      // let g1beg = Ellipse::new(grad1_p0, ( 2.,15.+5.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g1beg,);
      // let g1end = Ellipse::new(grad1_p1, ( 2.,15.+0.),  33.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g1end,);
      // let g2beg = Ellipse::new(grad2_p0, ( 2.,15.+5.),  99.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_beg, None, &g2beg,);
      // let g2end = Ellipse::new(grad2_p1, ( 2.,15.+0.),  99.0_f64.to_radians());scene.stroke(&pstr,Affine::IDENTITY, &col_end, None, &g2end,);
    }
  }

  pub fn ddd<I>(scene:&mut Scene, center:impl Into<Point>,r0:f64,  arc_beg:f64, arc_len_deg:f64,
    jn:JoinWhere, delta_transit:f64,
    col_beg:Color,col_end:Color,
    w1:f64,w2:f64, dpi:f64,
    precision_deg_per_step:f64,
    dash_off_deg:f64,dash_iter_deg_o:Option<I>,   dbg:u8,
    ) where I:IntoIterator, I::Item:Borrow<f64> {
      let mut is_vis_draw  = false; // whether a visible dash is drawing to clamp its first partial draw to the next. Switches to off when an invisible dash is "drawing"
      let mut dash_partial = 0.; // use as dash offset for the next segment to hide the partially drawn part
      let mut carry_over:f64 = 0.; // if Δstep covers 2 dash segments, the 1st one will store the remainer it didn't cover here for the 2nd to pick it up

      let delta_deg   	= arc_len_deg *       delta_transit ; let delta_rad    = delta_deg   .to_radians();
      let skip_beg_deg	= arc_len_deg * (1. - delta_transit); let skip_beg_rad = skip_beg_deg.to_radians();

      // Line width
      let wavg = (w1 + w2) / 2.; let wavgpx = (wavg * dpi).round() / dpi;
      let w_delta_avg = (w2 - w1).abs() / 2.;
      let w1px = (w1 * dpi).round() / dpi;
      let w2px = (w2 * dpi).round() / dpi;
      // ((w1 + sign1 * w_step * r) * dpi).round() / dpi; //transition is gradual, so can't be pixel-precise to avoid steps

      let precision_rad_per_step = precision_deg_per_step.to_radians(); //0.5° = 0.00873
      let steps_f           	= arc_len_deg / precision_deg_per_step; //180° → 360 steps at 0.5° precision
      let steps_i           	= steps_f as i32;
      let steps_delta_f     	= delta_deg / precision_deg_per_step; // 180°*33% / 0.5 = 118.8
      let steps_delta_i     	= steps_delta_f as i32; let steps_delta_if = f64::from(steps_delta_i); // 118 whole steps for later iteration and drawing by small step
      let delta_covered_deg 	=                  steps_delta_if  * precision_deg_per_step;
      // let steps_delta_rem	=  steps_delta_f - steps_delta_if; //0.8 steps not covered by whole
      let delta_rem_deg     	= delta_deg - delta_covered_deg; let delta_rem_rad = delta_rem_deg.to_radians();

      let w_step = w_delta_avg / steps_delta_f; //12/2/45 0.13 to reach average
      let w_per_step_i:f64	= w_delta_avg / f64::from(steps_delta_i); //to reach average
      let px1 = 1./r0; //length of 1px in rad
      let mut step_gap = if dbg>=1 {0.} else {px1}; //fix conflation artifacts outside of debug by overlapping segments (except the last one)

      let col_avg = col_beg.lerp(col_end,0.5,Default::default());

      let arc_len_rad = arc_len_deg.to_radians();
      let join_beg:f64 = match jn {
        JoinWhere::Beg	=> arc_beg,
        JoinWhere::End	=> arc_beg + (arc_len_deg - delta_deg).to_radians(),
      };

      let sign = match jn {
        JoinWhere::Beg	=> if w2px > wavg { 1.} else if w2px < wavg {-1.} else {0.}, //(from avg) ↑ if bigger, ↓ if smaller
        JoinWhere::End	=> if w1px > wavg {-1.} else if w1px < wavg { 1.} else {0.}, //(to   avg) ↓ if bigger, ↑ if smaller
      };
      let c = center.into(); let (cx,cy) = (c.x,c.y);

      // Dashes
      let (is_dash, dash_off_px, dash_off_rad, dash_off_pre_jn_rad, dash_iter_len_deg, dash_iter_len_rad, dash_iter_px, dash_iter_rad) = match dash_iter_deg_o {
        Some(dash_iter_deg)	=> {
          let deg_len = f64c::TAU * r0 / 360.; //2π*100/360 = 1.74
          let dash_off_rad     	= dash_off_deg.to_radians();
          let dash_off_px      	= dash_off_deg * deg_len;
          let dash_iter_px     	= dash_iter_deg.into_iter().map(|w|f_round(*w.borrow() * deg_len,6)).collect::<Vec<f64>>();
          let dash_iter_len_px 	= dash_iter_px.iter().sum::<f64>(); //5
          let dash_iter_rad    	= dash_iter_px.iter().map(|w| w / r0).collect::<Vec<f64>>();
          let dash_iter_len_deg	= dash_iter_len_px / deg_len; let dash_iter_len_rad = dash_iter_len_deg.to_radians(); //2.87356° 0.05015
          let dash_off_pre_jn_rad = match jn { // if join is after dashed line, that line can be cut and thus dash continuation should start from an adjusted offset
            JoinWhere::Beg	=> 0.,
            JoinWhere::End	=> skip_beg_rad % dash_iter_len_rad,
          };
          (true, dash_off_px, dash_off_rad, dash_off_pre_jn_rad, dash_iter_len_deg, dash_iter_len_rad, dash_iter_px, dash_iter_rad)
        },
        None	=> (false, 0.,0.,0., precision_deg_per_step,precision_rad_per_step, vec![precision_deg_per_step],vec![precision_rad_per_step]),
      };

      let dbg_col_beg = if dbg>=1{css::DARK_RED  }else{col_end};
      let dbg_col_end = if dbg>=1{css::DARK_GREEN}else{col_beg};
      let (grad_p0,grad_p1, col_stops) = match jn {
        JoinWhere::Beg	=> {
          let grad_beg_p0 = ( cx + r0*f64::cos(arc_beg                           ) , cy + r0*f64::sin(arc_beg                         ) );
          let grad_beg_p1 = ( cx + r0*f64::cos(arc_beg + delta_deg.to_radians()  ) , cy + r0*f64::sin(arc_beg + delta_deg.to_radians()) );
          let col_stops_beg = [
            (0.                                           ,col_avg    ),
            ((delta_deg  .to_radians()/arc_len_rad) as f32,col_end    ),
            ((delta_deg  .to_radians()/arc_len_rad) as f32,dbg_col_beg),
            ((arc_len_deg.to_radians()/arc_len_rad) as f32,dbg_col_beg),
            ];
          (grad_beg_p0,grad_beg_p1, col_stops_beg)},
        JoinWhere::End	=> {
          let grad_end_p0 = ( cx + r0*f64::cos(arc_beg + skip_beg_rad            ) , cy + r0*f64::sin(arc_beg + skip_beg_rad          ));
          let grad_end_p1 = ( cx + r0*f64::cos(arc_beg + arc_len_rad             ) , cy + r0*f64::sin(arc_beg + arc_len_rad           ));
          let col_stops_end = [
            (0.                                           ,dbg_col_end),
            ((skip_beg_rad            /arc_len_rad) as f32,dbg_col_end),
            ((skip_beg_rad            /arc_len_rad) as f32,col_beg    ),
            ((arc_len_deg.to_radians()/arc_len_rad) as f32,col_avg    ),
            ];
          (grad_end_p0,grad_end_p1, col_stops_end)},
      };
      let grad = Gradient::new_sweep((cx,cy),arc_beg as f32,(arc_beg + arc_len_rad) as f32).with_stops(col_stops);

      if let JoinWhere::End = jn {// Segment 1: ~join part is 2nd (at the end)
        // Draw pre-gradwidth segment separately without the extra iterator
        let c = Arc::new((cx,cy), (r0,r0)   ,  arc_beg,skip_beg_rad, 0.);
        let stroke_c = get_stroke_end(w1px);
        scene.stroke(&stroke_c, Affine::IDENTITY, &dbg_col_end, None, &c,);
      }

      // TODO: force precision to be so that one step is never bigger than a dash set length, otherwise would need to repeat the full dash parsing logic here
        // or maybe have a better loop logic that can handle it?
      let is_extra_step = delta_rem_rad > 0.;
      let steps_delta_xt = if is_extra_step {steps_delta_i} else {steps_delta_i-1};
      let w_beg = match jn {
        JoinWhere::Beg	=> wavg,
        JoinWhere::End	=> w1px,};
      let w_end = match jn {
        JoinWhere::Beg	=> w2px,
        JoinWhere::End	=> wavg,};

      // Calculate which of the steps will match the last visible dash, and then do NOT extend the last step to cover for occlusion artifacts of other steps
      let (step_ix, dash_vis_ix) = if !is_dash {
        let arc_in_dashes = (arc_len_rad + dash_off_rad) % dash_iter_len_rad; // last partial arc's space that would hold a dash (in arc's coordinates)
        let last_dash_len = if arc_in_dashes > 0. {arc_in_dashes} else {dash_iter_len_rad}; // last full or partial arc's space that would hold a dash
        let mut dash_ix 	= 0;
        let mut dash_rem	= dash_iter_rad[dash_ix];
        let mut is_active = true;
        while dash_rem < last_dash_len { // Find place in dashes array that covers arc's end
          dash_ix = (dash_ix + 1) % dash_iter_rad.len();
          dash_rem += dash_iter_rad[dash_ix];
          is_active = !is_active;
        }
        let dash_vis_ix	= dash_ix - if is_active {0} else {1}; // last visible dash in at the end of the line
        let mut dash_last_vis_end = dash_iter_rad[dash_vis_ix]; // end of last vis dash = length of all dashes up to and including the last visible one
        if dash_vis_ix > 0 {for di in 1..dash_vis_ix {dash_last_vis_end += dash_iter_rad[di-1]}};
        let space_leftover = (last_dash_len - dash_last_vis_end).max(0.);
        let dash_last_vis_end_arc = arc_len_rad - space_leftover; // same, but in arc's coordinates
        let space_leftover_arc = dash_last_vis_end_arc - dash_last_vis_end_arc;

        // Get the last step that covers the last visible dash
        let step_ix = if is_extra_step && space_leftover_arc < delta_rem_rad {f64::from(steps_delta_xt)} else {
          dash_last_vis_end_arc.div_euclid(precision_rad_per_step) + if dash_last_vis_end_arc % precision_rad_per_step > 0. {1.}else{0.}
        };
        if dbg>=2{let _six1 = step_ix + 1.0; let _steps=steps_delta_i; let _dix1 = dash_vis_ix + 1; let _dashes=dash_iter_px.len();
          println!("step_ix={_six1}¦{_steps}{} dash_vis_ix={_dix1}¦{_dashes} arc_len={}° ╍off={} ╍len={dash_iter_len_deg}° arc%╍={}°=({}+{})%{} ␠leftover={}°"
          ,if is_extra_step {"+1"}else{""}
          ,arc_len_rad.to_degrees(), dash_off_rad.to_degrees(), arc_in_dashes.to_degrees(),dash_off_rad.to_degrees(),arc_len_rad.to_degrees(),dash_iter_len_rad.to_degrees()
          ,space_leftover.to_degrees());}
        (step_ix as usize, dash_vis_ix)
      } else {(0usize, 0usize)};


      for i in 0..=steps_delta_xt { let r = f64::from(i); let is_last = i == steps_delta_xt;
        // NB! last step needs special handling since it's a fractional one, so not full "precision length"!
        let step_width = if is_last && is_extra_step	{delta_rem_rad
        } else                                      	{precision_rad_per_step};
        let seg0 = r * precision_rad_per_step; // segment beg in our arc coords (arc start = 0), last regular starts at the same old width, but itself will have a smaller step_width
        let rad0 = join_beg + seg0;
        let rad1 = rad0 + step_width; // todo debug only
        // let c = Arc::new((cx,cy), (r0,r0) ,  rad0,step_width+gap_correct, 0.); //arc bugs with gaps
        // let c = CircleSegment::new((cx,cy), r0,r0   ,  rad0,step_width); // alt fix
        //                          center  rout/in    ∠start ∠sweep
        let cw = if is_last	{w_end
        } else             	{w_beg + sign * r * w_per_step_i};

        let stroke_c = if dbg>=2 {match jn {
          JoinWhere::Beg	=> get_stroke_end(w2px),
          JoinWhere::End	=> get_stroke_end(w1px),}
        } else{get_stroke_end(cw)};

        let seg_off =  dash_off_pre_jn_rad + dash_off_rad         % dash_iter_len_rad;
        let seg_be_ = (dash_off_pre_jn_rad + dash_off_rad + seg0) % dash_iter_len_rad;// cut off arc that fit into the previous dash set, so this is our segment beginning in the coordinate system of a dash set (set's begin = 0)
        let (seg_beg,seg_count) = if (dash_iter_len_rad - seg_be_).abs() < epsi { // segment starts at dash set end's, so shift it to next dash set's begin (float imprecision prevents clean division)
          (0.     , (dash_off_pre_jn_rad + dash_off_rad + seg0).div_euclid(dash_iter_len_rad) + 1.+1.)
        } else {
          (seg_be_, (dash_off_pre_jn_rad + dash_off_rad + seg0).div_euclid(dash_iter_len_rad) + 1.   )
        };
        let seg_count = (dash_off_pre_jn_rad + dash_off_rad + seg0).div_euclid(dash_iter_len_rad) + 1.;
        let seg_end = seg_beg + step_width;
        let mut d_beg = 0.; // length up to the beginning of this dash = ∑ of all previous dash lens

        //     ──────  —————  outer line dash pattern (todo: what if our line is bigger than 1 pattern? like here)
        //        ┌─────┐     our   line (always starts later due to "rad0 % dash_iter_len_rad")
        // seg_beg┘     └seg_end
        //         ↑↑  ↑ draw, overlaps with   active
        //           ↑↑  skip, overlaps with inactive
        if dbg>=2 {if i == 0 {println!("\n\n—————Σⁱ={steps_delta_xt: >3}——╍╍№{} Σ{dash_iter_len_deg: >4.1}° off{dash_off_deg: >4.1}° {dash_iter_px:?}°╍╍——beg {: >4.1}° Δ{delta_covered_deg: >4.1}° + {delta_rem_deg: >4.1}° rem = Δ{delta_deg: >4.1}°————————————————————————————————————————————"
          ,dash_iter_px.len(), arc_beg.to_degrees())}};
        let mut dr = 0; // track dash 🗘
        let mut step_covered = step_width; // track Σ dash_iter_len_rad covering each Δstep
        while step_covered > 0.  {
          dr += 1;
          step_covered -= dash_iter_len_rad;
        let mut j = 0;
        let mut is_visible = false;
        let mut prev_draw_len:f64 = 0.; // store a sum of previously drawn dashes (vis+invis) so that if this dash doesn't cover the full Δstep or Σdash_len, we can see which part of it was covered before and which should go as Δover to the next step
        // if seg_count == 1. {

        let mut dash_ix: usize = 0;
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
              if dbg>=1	{scene.stroke(&stroke_c, Affine::IDENTITY, css::MAGENTA , None, &c,);
              } else   	{scene.stroke(&stroke_c, Affine::IDENTITY, &grad        , None, &c,);}
              carry_over = 0.;
              dbgprint = true;
              if dbg>=4 && dbgprint {
              println!("{i} drawing Δover {: >2.1} @ {: >3.2} = ({: >2.1}-{: >2.1}-Δ{: >2.1}) → {: >3.2}",carry_over.to_degrees()
                ,(rad1 - step_width - carry_over).to_degrees(),rad1.to_degrees(),step_width.to_degrees(),carry_over.to_degrees()
                ,(rad1 - step_width).to_degrees());}
            }
            if draw_len > 0.0 { // 1st draw starts @ seg end to attach to the next draw in case of partials
              prev_draw_len += draw_len;
              if dbg == 0 && step_gap > 0. && (i as usize) == step_ix && dash_ix == dash_vis_ix {step_gap=0.}; // don't bleed the last dash's visible end into the next segment. Deal with gap between arcs by drawing the next arc earlier?
              let c = if is_vis_draw   {Arc::new((cx,cy), (r0,r0)   ,rad0         ,draw_len + step_gap, 0.)
              } else {is_vis_draw=true; Arc::new((cx,cy), (r0,r0)   ,rad1-draw_len,draw_len + step_gap, 0.)};
              if dbg>=1 {
                if is_last	{scene.stroke(&stroke_c, Affine::IDENTITY, &css::LIME, None, &c,);
                } else    	{scene.stroke(&stroke_c, Affine::IDENTITY, &grad     , None, &c,);}
              } else      	{scene.stroke(&stroke_c, Affine::IDENTITY, &grad     , None, &c,);}
              if is_last && draw_len < *dash_i - epsi { // drawn something, but not the full visible dash
                let part_len = draw_end - d_beg; //how much of an existing dash is covered by all draws, incl. last
                dash_partial = (d_beg + part_len) * r0; // add all prior dash segments within a set
                // println!("! last +visible +draw ╍w {: >.4}° − {: >.4}° actual = {: >.4}° {: >.4}°part_len partial ({:.1}px) rad1 {:.3}°"
                //   ,dash_i.to_degrees(),draw_len.to_degrees(), (d_end.min(seg_end) - d_beg).to_degrees(),part_len.to_degrees(), dash_partial,rad1.to_degrees());
              }
            } else {is_vis_draw=false;}
            // if rad0       <=       d_end
            // &&    seg_end >= d_beg  { // (alt check) our segment overlaps with this dash
            if dbg>=2 && (dbgprint || i == 0 || is_last || (78<= i && i <=83)) {println!( //👁👀👓
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
                  if dbg>=1	{scene.stroke(&stroke_c, Affine::IDENTITY, css::MAGENTA , None, &c,);
                  } else   	{scene.stroke(&stroke_c, Affine::IDENTITY, &grad        , None, &c,);}
                  dash_partial = carry_over_delta * r0; carry_over = 0.;
                  // println!("last step - drawn next dash since it won't be handled later!");
                // } else {println!("  Δover {: >4.1}° = step_w {: >4.1}° - {: >4.1}° drawn",carry_over.to_degrees(),step_width.to_degrees(),draw_len.to_degrees());
                }
              }
            }
            if is_last { // last invisible dash also needs to signal its width to update offset of the next arc
              let part_len = draw_end - d_beg; //how much of an existing dash is covered by all draws, incl. last
              if   draw_len > 0. //drawn something… ↙some float rounding error
                && part_len < *dash_i - 0.00000000001 { //…but not the full invisible dash
                dash_partial = (d_beg + part_len) * r0; //≝draw_end add all prior dash segments within a set
                // println!("{}№{seg_count} last -visible +draw ╍beg {: >3.3}° draw_end {: >3.3}°≝partial ({:.1}px) Δ{: >3.3}° Δstep {: >3.3}° drawn │ ╍w {: >.2}° left {: >.2}° rad1 {:.3}°"
                //   ,if dash_partial > 0. {"✓"}else{"✗"}
                //   ,d_beg.to_degrees(),draw_end.to_degrees(),dash_partial,part_len.to_degrees(),draw_len.to_degrees()
                //   ,dash_i.to_degrees(),(dash_i-part_len).to_degrees(),rad1.to_degrees());
              }
            }
            if dbg>=3 && (dbgprint || i == 0 || is_last || (78<= i && i <=83)) {println!( //👁👀👓
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
          dash_ix += 1;
        }
        }
      }

      if let JoinWhere::Beg = jn { // Segment 2: ~join part is 1st (at the start)
        // Draws pos-gradwidth segment separately without the extra iterator, including leftovers from whole steps not covering the full range
        // println!("2nd curve@end: {: >3.0}° + Δ{: >3.0}° ⇒ {: >3.0}° + {: >3.0} skip_beg ╍part={: >3.0}° ({dash_partial: >3.0}px)"
        //   ,arc_beg.to_degrees(),delta_rad.to_degrees(),(arc_beg + delta_rad).to_degrees(),skip_beg_rad.to_degrees(),(dash_partial/r0).to_degrees());
        let c = Arc::new((cx,cy), (r0,r0) , arc_beg + delta_rad, skip_beg_rad, 0.);
        // let stroke_c = get_stroke_end(w2px); // todo↓ make dashes optional
        let stroke_c = if is_dash { // use remainder from the previous segment so that the total matches the style as though it were drawn in one step
          get_stroke_end      (w2px).with_dashes(dash_partial,&dash_iter_px)
        } else {get_stroke_end(w2px)};
        if dbg>=1	{scene.stroke(&stroke_c, Affine::IDENTITY, &css::DARK_RED, None, &c,);
        } else   	{scene.stroke(&stroke_c, Affine::IDENTITY, &col_end      , None, &c,);}
      }

    if dbg >=1 { // DEBUG copy smaller (including dashes, should perfectly align as dash length/offsets are adjusted per difference in size)
      let r00 = match jn {
        JoinWhere::Beg	=> r0 - (w2px + w2px.max(wavgpx)) / 2.,
        JoinWhere::End	=> r0 - (w1px + w1px.max(wavgpx)) / 2.,};
      let wpx = match jn {
        JoinWhere::Beg	=> w2px,
        JoinWhere::End	=> w1px,};
      let grad = Gradient::new_sweep((cx,cy),             arc_beg as f32,(arc_beg + arc_len_rad) as f32).with_stops(col_stops);
      let c    = Arc::new           ((cx,cy), (r00,r00),  arc_beg       ,           arc_len_rad, 0.);
      let stroke_c = if is_dash {
        get_stroke_end(wpx).with_dashes(dash_off_px*r00/r0,dash_iter_px.iter().map(|w|f_round(w*r00/r0,6)).collect::<Vec<f64>>())
      } else {get_stroke_end(wpx)};
      scene.stroke(&stroke_c, Affine::IDENTITY, &grad, None, &c,);
    }
  }
}