function [fine_model,fine_p, fine_sW] = interpolate(fine_model,coarse_model, coarse_p_ad, coarse_sW_ad)
%%Apply the coarse grid partition field to interpolate the variables from the
%%coarse grid into the fine grid

  %fine_p_ad = coarseDataToFine(CG, coarse_p_ad);
  
  %fine_sW_ad = coarseDataToFine(CG, coarse_sW_ad);
  
  fine_p = coarse_p_ad.val(coarse_model.grid.partition);
  
  fine_sW = coarse_sW_ad.val(coarse_model.grid.partition);
  
  fine_model.cycle = coarse_model.cycle;
end