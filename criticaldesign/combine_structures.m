%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   This file is part of dAEDalus criticaldesign
%                   Copyright (C) 2011, Klaus Seywald
%     Author:   	Klaus Seywald
%                   klaus.seywald@mytum.de
%                   seywald@kth.se

function  [wingstructure] = combine_structures(wingstructure_prv,wingstructure,index)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 if isempty(wingstructure.beam)

    for i=1:length(wingstructure.beamelement)
        
        if isa(wingstructure.beamelement(i).crosssection,'class_crosssection_fuselage')
            if(wingstructure_prv.beamelement(i).crosssection.t_sk_eq>=wingstructure.beamelement(i).crosssection.t_sk_eq);
                wingstructure.beamelement(i).crosssection.t_sk_eq=wingstructure_prv.beamelement(i).crosssection.t_sk_eq;
                wingstructure.beamelement(i).crosssection.t_sk_eq_lc_idx=wingstructure_prv.beamelement(i).crosssection.t_sk_eq_lc_idx;
            else
               wingstructure.beamelement(i).crosssection.t_sk_eq_lc_idx=index; 
            end
        
            if(wingstructure_prv.beamelement(i).crosssection.t_fr_eq>=wingstructure.beamelement(i).crosssection.t_fr_eq);
                wingstructure.beamelement(i).crosssection.t_fr_eq=wingstructure_prv.beamelement(i).crosssection.t_fr_eq;
                wingstructure.beamelement(i).crosssection.t_fr_eq_lc_idx=wingstructure_prv.beamelement(i).crosssection.t_fr_eq_lc_idx;
            else
               wingstructure.beamelement(i).crosssection.t_fr_eq_lc_idx=index; 
            end
            
        else
        
            if(wingstructure_prv.beamelement(i).crosssection.t_sp_fr>=wingstructure.beamelement(i).crosssection.t_sp_fr);
                wingstructure.beamelement(i).crosssection.t_sp_fr=wingstructure_prv.beamelement(i).crosssection.t_sp_fr;
                wingstructure.beamelement(i).crosssection.t_sp_fr_lc_idx=wingstructure_prv.beamelement(i).crosssection.t_sp_fr_lc_idx;
            else
               wingstructure.beamelement(i).crosssection.t_sp_fr_lc_idx=index; 
            end

            if(wingstructure_prv.beamelement(i).crosssection.t_sp_re>=wingstructure.beamelement(i).crosssection.t_sp_re);
                wingstructure.beamelement(i).crosssection.t_sp_re=wingstructure_prv.beamelement(i).crosssection.t_sp_re;
                wingstructure.beamelement(i).crosssection.t_sp_re_lc_idx=wingstructure_prv.beamelement(i).crosssection.t_sp_re_lc_idx;
            else
               wingstructure.beamelement(i).crosssection.t_sp_re_lc_idx=index; 
            end

            if(wingstructure_prv.beamelement(i).crosssection.t_sk_up>=wingstructure.beamelement(i).crosssection.t_sk_up);
                wingstructure.beamelement(i).crosssection.t_sk_up=wingstructure_prv.beamelement(i).crosssection.t_sk_up;
                wingstructure.beamelement(i).crosssection.t_sk_up_lc_idx=wingstructure_prv.beamelement(i).crosssection.t_sk_up_lc_idx;
            else
               wingstructure.beamelement(i).crosssection.t_sk_up_lc_idx=index; 
            end  

            if(wingstructure_prv.beamelement(i).crosssection.t_sk_lo>=wingstructure.beamelement(i).crosssection.t_sk_lo);
                wingstructure.beamelement(i).crosssection.t_sk_lo=wingstructure_prv.beamelement(i).crosssection.t_sk_lo;
                wingstructure.beamelement(i).crosssection.t_sk_lo_lc_idx=wingstructure_prv.beamelement(i).crosssection.t_sk_lo_lc_idx;
            else
               wingstructure.beamelement(i).crosssection.t_sk_lo_lc_idx=index; 
            end
        end
        %recalculate crosssectional parameters
        wingstructure.beamelement(i)=wingstructure.beamelement(i).f_calcCrossProp();
    end
    wingstructure.update_K=1;
 else
   for j=1:length(wingstructure.beam)
     for i=1:length(wingstructure.beam(j).beamelement)
         
        if isa(wingstructure.beam(j).beamelement(i).crosssection,'class_crosssection_fuselage')
            if(wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_eq>=wingstructure.beam(j).beamelement(i).crosssection.t_sk_eq);
                wingstructure.beam(j).beamelement(i).crosssection.t_sk_eq=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_eq;
                wingstructure.beam(j).beamelement(i).crosssection.t_sk_eq_lc_idx=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_eq_lc_idx;
            else
               wingstructure.beam(j).beamelement(i).crosssection.t_sk_eq_lc_idx=index; 
            end

            if(wingstructure_prv.beam(j).beamelement(i).crosssection.t_fr_eq>=wingstructure.beam(j).beamelement(i).crosssection.t_fr_eq);
                wingstructure.beam(j).beamelement(i).crosssection.t_fr_eq=wingstructure_prv.beam(j).beamelement(i).crosssection.t_fr_eq;
                wingstructure.beam(j).beamelement(i).crosssection.t_fr_eq_lc_idx=wingstructure_prv.beam(j).beamelement(i).crosssection.t_fr_eq_lc_idx;
            else
               wingstructure.beam(j).beamelement(i).crosssection.t_fr_eq_lc_idx=index; 
            end
        elseif isa(wingstructure.beam(j).beamelement(i).crosssection,'class_crosssection_wingbox')            
        
            if(wingstructure_prv.beam(j).beamelement(i).crosssection.t_sp_fr>=wingstructure.beam(j).beamelement(i).crosssection.t_sp_fr);
                wingstructure.beam(j).beamelement(i).crosssection.t_sp_fr=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sp_fr;
                wingstructure.beam(j).beamelement(i).crosssection.t_sp_fr_lc_idx=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sp_fr_lc_idx;
            else
               wingstructure.beam(j).beamelement(i).crosssection.t_sp_fr_lc_idx=index; 
            end

            if(wingstructure_prv.beam(j).beamelement(i).crosssection.t_sp_re>=wingstructure.beam(j).beamelement(i).crosssection.t_sp_re);
                wingstructure.beam(j).beamelement(i).crosssection.t_sp_re=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sp_re;
                wingstructure.beam(j).beamelement(i).crosssection.t_sp_re_lc_idx=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sp_re_lc_idx;
            else
               wingstructure.beam(j).beamelement(i).crosssection.t_sp_re_lc_idx=index; 
            end

            if(wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_up>=wingstructure.beam(j).beamelement(i).crosssection.t_sk_up);
                wingstructure.beam(j).beamelement(i).crosssection.t_sk_up=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_up;
                wingstructure.beam(j).beamelement(i).crosssection.t_sk_up_lc_idx=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_up_lc_idx;
            else
               wingstructure.beam(j).beamelement(i).crosssection.t_sk_up_lc_idx=index; 
            end  

            if(wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_lo>=wingstructure.beam(j).beamelement(i).crosssection.t_sk_lo);
                wingstructure.beam(j).beamelement(i).crosssection.t_sk_lo=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_lo;
                wingstructure.beam(j).beamelement(i).crosssection.t_sk_lo_lc_idx=wingstructure_prv.beam(j).beamelement(i).crosssection.t_sk_lo_lc_idx;
            else
               wingstructure.beam(j).beamelement(i).crosssection.t_sk_lo_lc_idx=index; 
            end
        end
        %recalculate crosssectional parameters
        wingstructure.beam(j).beamelement(i)=wingstructure.beam(j).beamelement(i).f_calcCrossProp();
    end
    wingstructure.beam(j).update_K=1;
    wingstructure.beam(j).update_M=1;
  end
end
        
end