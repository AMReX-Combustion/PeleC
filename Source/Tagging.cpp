#include <AMReX_ParmParse.H>

#include "PeleC.H"
#include "Tagging.H"

void
PeleC::read_tagging_params()
{
  const std::string tag_prefix = "tagging";
  amrex::ParmParse pp(tag_prefix);

  pp.query("denerr", tagging_parm->denerr);
  pp.query("max_denerr_lev", tagging_parm->max_denerr_lev);
  pp.query("dengrad", tagging_parm->dengrad);
  pp.query("max_dengrad_lev", tagging_parm->max_dengrad_lev);
  pp.query("denratio", tagging_parm->denratio);
  pp.query("max_denratio_lev", tagging_parm->max_denratio_lev);

  pp.query("presserr", tagging_parm->presserr);
  pp.query("max_presserr_lev", tagging_parm->max_presserr_lev);
  pp.query("pressgrad", tagging_parm->pressgrad);
  pp.query("max_pressgrad_lev", tagging_parm->max_pressgrad_lev);

  pp.query("velerr", tagging_parm->velerr);
  pp.query("max_velerr_lev", tagging_parm->max_velerr_lev);
  pp.query("velgrad", tagging_parm->velgrad);
  pp.query("max_velgrad_lev", tagging_parm->max_velgrad_lev);

  pp.query("vorterr", tagging_parm->vorterr);
  pp.query("max_vorterr_lev", tagging_parm->max_vorterr_lev);

  pp.query("temperr", tagging_parm->temperr);
  pp.query("max_temperr_lev", tagging_parm->max_temperr_lev);
  pp.query("lotemperr", tagging_parm->lotemperr);
  pp.query("max_lotemperr_lev", tagging_parm->max_lotemperr_lev);
  pp.query("tempgrad", tagging_parm->tempgrad);
  pp.query("max_tempgrad_lev", tagging_parm->max_tempgrad_lev);

  pp.query("ftracerr", tagging_parm->ftracerr);
  pp.query("max_ftracerr_lev", tagging_parm->max_ftracerr_lev);
  pp.query("ftracgrad", tagging_parm->ftracgrad);
  pp.query("max_ftracgrad_lev", tagging_parm->max_ftracgrad_lev);

  pp.query("vfracerr", tagging_parm->vfracerr);
  pp.query("max_vfracerr_lev", tagging_parm->max_vfracerr_lev);

  pp.query("eb_refine_type", tagging_parm->eb_refine_type);
  pp.query("max_eb_refine_lev", tagging_parm->max_eb_refine_lev);
  pp.query("eb_detag_factor", tagging_parm->detag_eb_factor);
  tagging_parm->detag_eb_factor *= std::sqrt(2.0);
  if (
    (tagging_parm->eb_refine_type != "static") &&
    (tagging_parm->eb_refine_type != "adaptive")) {
    amrex::Abort("refine_eb_type can only be 'static' or 'adaptive'");
  }
  if (tagging_parm->eb_refine_type == "adaptive") {
    tagging_parm->min_eb_refine_lev = 0;
    pp.query("min_eb_refine_lev", tagging_parm->min_eb_refine_lev);
    tagging_parm->adapt_eb_refined_lev = tagging_parm->min_eb_refine_lev;
  }

  // amrex tagging
  amrex::Vector<std::string> refinement_indicators;
  pp.queryarr(
    "refinement_indicators", refinement_indicators, 0,
    pp.countval("refinement_indicators"));
  for (const auto& refinement_indicator : refinement_indicators) {
    std::string ref_prefix = tag_prefix + "." + refinement_indicator;
    amrex::ParmParse ppr(ref_prefix);

    // Tag a given box
    amrex::RealBox realbox;
    if (ppr.countval("in_box_lo") > 0) {
      amrex::Vector<amrex::Real> box_lo(AMREX_SPACEDIM);
      amrex::Vector<amrex::Real> box_hi(AMREX_SPACEDIM);
      ppr.getarr("in_box_lo", box_lo, 0, static_cast<int>(box_lo.size()));
      ppr.getarr("in_box_hi", box_hi, 0, static_cast<int>(box_hi.size()));
      realbox = amrex::RealBox(box_lo.data(), box_hi.data());
    }

    amrex::AMRErrorTagInfo info;

    if (realbox.ok()) {
      info.SetRealBox(realbox);
    }

    if (ppr.countval("start_time") > 0) {
      amrex::Real min_time;
      ppr.get("start_time", min_time);
      info.SetMinTime(min_time);
    }

    if (ppr.countval("end_time") > 0) {
      amrex::Real max_time;
      ppr.get("end_time", max_time);
      info.SetMaxTime(max_time);
    }

    if (ppr.countval("max_level") > 0) {
      int tag_max_level;
      ppr.get("max_level", tag_max_level);
      info.SetMaxLevel(tag_max_level);
    }

    bool itexists = false;
    int index = State_Type;
    int scomp = 0;
    if (ppr.countval("value_greater") > 0) {
      amrex::Real value;
      ppr.get("value_greater", value);
      std::string field;
      ppr.get("field_name", field);
      tagging_parm->err_tags.push_back(
        amrex::AMRErrorTag(value, amrex::AMRErrorTag::GREATER, field, info));
      itexists =
        derive_lst.canDerive(field) || isStateVariable(field, index, scomp);
    } else if (ppr.countval("value_less") > 0) {
      amrex::Real value;
      ppr.get("value_less", value);
      std::string field;
      ppr.get("field_name", field);
      tagging_parm->err_tags.push_back(
        amrex::AMRErrorTag(value, amrex::AMRErrorTag::LESS, field, info));
      itexists =
        derive_lst.canDerive(field) || isStateVariable(field, index, scomp);
    } else if (ppr.countval("vorticity_greater") > 0) {
      amrex::Real value;
      ppr.get("vorticity_greater", value);
      const std::string field = "magvort";
      tagging_parm->err_tags.push_back(
        amrex::AMRErrorTag(value, amrex::AMRErrorTag::VORT, field, info));
      itexists =
        derive_lst.canDerive(field) || isStateVariable(field, index, scomp);
    } else if (ppr.countval("adjacent_difference_greater") > 0) {
      amrex::Real value;
      ppr.get("adjacent_difference_greater", value);
      std::string field;
      ppr.get("field_name", field);
      tagging_parm->err_tags.push_back(
        amrex::AMRErrorTag(value, amrex::AMRErrorTag::GRAD, field, info));
      itexists =
        derive_lst.canDerive(field) || isStateVariable(field, index, scomp);
    } else if (realbox.ok()) {
      tagging_parm->err_tags.push_back(amrex::AMRErrorTag(info));
      itexists = true;
    } else {
      amrex::Abort(
        "Unrecognized refinement indicator for " + refinement_indicator);
    }

    if (!itexists) {
      amrex::Error(
        "PeleC::read_tagging_params(): unknown variable field for criteria " +
        refinement_indicator);
    }
  }
}
