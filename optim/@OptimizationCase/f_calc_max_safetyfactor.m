function [obj,varargout] = f_calc_max_safetyfactor(obj,nomDisp,nomLoads,beamID,pos,test_maxgusts)

nT=size(obj.aeSSM.gustData.modes,1); %number of time step
nL=size(obj.aeSSM.gustData.modes,3); %number of gust lengths
nEl=length(obj.str.beam(beamID).beamelement); 

%matrix with computed values of failure index 
allFailIdx_fs = zeros(nT,nEl,nL);
allFailIdx_rs = zeros(nT,nEl,nL);
allFailIdx_us = zeros(nT,nEl,nL);
allFailIdx_ls = zeros(nT,nEl,nL);


% displacement are computed depending on the critical case (positive or
% negative)
disp = zeros(nT,size(obj.aeSSM.strModel.modalBase,1),nL);
if pos
    for iT=1:nT
        for iL=1:nL
            disp(iT,:,iL)=nomDisp+obj.aeSSM.strModel.modalBase*obj.aeSSM.gustData.modes(iT,:,iL)';
        end
    end
else
    for iT=1:nT
        for iL=1:nL
            disp(iT,:,iL)=nomDisp-obj.aeSSM.strModel.modalBase*obj.aeSSM.gustData.modes(iT,:,iL)';
        end
    end
end

% remove the calculation of the failure index so only strains are computed
obj.str.beam(beamID).calcul_failIdx = 0;

jT = 10; %To change to a value dependent on nT & nL ?
jGl = 5;

sample_time = 1:jT:nT;
sample_gl = 1:jGl:nL;

IdxEstimation_fs = zeros(nEl,nT,nL);
IdxEstimation_rs = zeros(nEl,nT,nL);
IdxEstimation_us = zeros(nEl,nT,nL);
IdxEstimation_ls = zeros(nEl,nT,nL);


% estimation on a scarce grid of the failure index
tic;
for iT=sample_time
    for iL=sample_gl

        obj.str.nodal_deflections=disp(iT,:,iL)';
        obj.str=obj.str.f_postprocess();
        obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
        
        for iBEl = 1:length(obj.str.beam(beamID).beamelement)
                      
            u1 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_fs.u1;
            u2 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_fs.u2;
            u3 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_fs.u3;
            u4 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_fs.u4;
            u5 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_fs.u5;
            u6 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_fs.u6;
            
            i1 = 0; i2 = sqrt((u4^2+2*u1+abs(u4)*sqrt(u4^2+4*u1))/(2*u1^2));
            tst = 4*u6^2*i2.^2 - 4*u6*u1*i2.^2 + 4*(1 - u2*i1 - u3*i1.^2)*(u1 - u6) + (u4 + u5*i1).^2;
            
            S = obj.str.beam(beamID).beamelement(iBEl).crosssection.strain_fs(1:3,:);
            I1 = S(1,:) + S(2,:);
            I2 = sqrt(((S(1,:)-S(2,:))/2).^2 + S(3,:).^2);
            
            if tst >= 0
                IdxEstimation_fs(iBEl,iT,iL) = max(4*u6^2*I2.^2 - 4*u6*u1*I2.^2 + 4*(1 - u2*I1 - u3*I1.^2)*(u1 - u6) + (u4 + u5*I1).^2);
            else
                IdxEstimation_fs(iBEl,iT,iL) = max(u1.^2*I2.^4 - I2.^2.*(u4 + u5*I1).^2 - 2*u1*I2.^2.*(1-u2*I1-u3*I1.^2)+(1-u2*I1-u3*I1.^2).^2);
            end
            
        end  
        for iBEl = 1:length(obj.str.beam(beamID).beamelement)
            u1 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_rs.u1;
            u2 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_rs.u2;
            u3 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_rs.u3;
            u4 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_rs.u4;
            u5 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_rs.u5;
            u6 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_rs.u6;
            
            i1 = 0; i2 = sqrt((u4^2+2*u1+abs(u4)*sqrt(u4^2+4*u1))/(2*u1^2));
            tst = 4*u6^2*i2.^2 - 4*u6*u1*i2.^2 + 4*(1 - u2*i1 - u3*i1.^2)*(u1 - u6) + (u4 + u5*i1).^2;
            
            S = obj.str.beam(beamID).beamelement(iBEl).crosssection.strain_rs(1:3,:);
            I1 = S(1,:) + S(2,:);
            I2 = sqrt(((S(1,:)-S(2,:))/2).^2 + S(3,:).^2);
            
            if tst >= 0
                IdxEstimation_rs(iBEl,iT,iL) = max(4*u6^2*I2.^2 - 4*u6*u1*I2.^2 + 4*(1 - u2*I1 - u3*I1.^2)*(u1 - u6) + (u4 + u5*I1).^2);
            else
                IdxEstimation_rs(iBEl,iT,iL) = max(u1.^2*I2.^4 - I2.^2.*(u4 + u5*I1).^2 - 2*u1*I2.^2.*(1-u2*I1-u3*I1.^2)+(1-u2*I1-u3*I1.^2).^2);
            end
            
            
        end
        for iBEl = 1:length(obj.str.beam(beamID).beamelement)
            u1 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_up.u1;
            u2 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_up.u2;
            u3 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_up.u3;
            u4 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_up.u4;
            u5 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_up.u5;
            u6 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_up.u6;
            
            i1 = 0; i2 = sqrt((u4^2+2*u1+abs(u4)*sqrt(u4^2+4*u1))/(2*u1^2));
            tst = 4*u6^2*i2.^2 - 4*u6*u1*i2.^2 + 4*(1 - u2*i1 - u3*i1.^2)*(u1 - u6) + (u4 + u5*i1).^2;
            
            S = obj.str.beam(beamID).beamelement(iBEl).crosssection.strain_sk_up(1:3,:);
            I1 = S(1,:) + S(2,:);
            I2 = sqrt(((S(1,:)-S(2,:))/2).^2 + S(3,:).^2);
            
            if tst >= 0
                IdxEstimation_us(iBEl,iT,iL) = max(4*u6^2*I2.^2 - 4*u6*u1*I2.^2 + 4*(1 - u2*I1 - u3*I1.^2)*(u1 - u6) + (u4 + u5*I1).^2);
            else
                IdxEstimation_us(iBEl,iT,iL) = max(u1.^2*I2.^4 - I2.^2.*(u4 + u5*I1).^2 - 2*u1*I2.^2.*(1-u2*I1-u3*I1.^2)+(1-u2*I1-u3*I1.^2).^2);
            end
            
            
        end
        for iBEl = 1:length(obj.str.beam(beamID).beamelement)
            u1 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_lo.u1;
            u2 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_lo.u2;
            u3 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_lo.u3;
            u4 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_lo.u4;
            u5 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_lo.u5;
            u6 = obj.str.beam(beamID).beamelement(iBEl).crosssection.laminate_sk_lo.u6;
            
            i1 = 0; i2 = sqrt((u4^2+2*u1+abs(u4)*sqrt(u4^2+4*u1))/(2*u1^2));
            tst = 4*u6^2*i2.^2 - 4*u6*u1*i2.^2 + 4*(1 - u2*i1 - u3*i1.^2)*(u1 - u6) + (u4 + u5*i1).^2;
            
            S = obj.str.beam(beamID).beamelement(iBEl).crosssection.strain_sk_lo(1:3,:);
            I1 = S(1,:) + S(2,:);
            I2 = sqrt(((S(1,:)-S(2,:))/2).^2 + S(3,:).^2);
            
            if tst >= 0
                IdxEstimation_ls(iBEl,iT,iL) = max(4*u6^2*I2.^2 - 4*u6*u1*I2.^2 + 4*(1 - u2*I1 - u3*I1.^2)*(u1 - u6) + (u4 + u5*I1).^2);
            else
                IdxEstimation_ls(iBEl,iT,iL) = max(u1.^2*I2.^4 - I2.^2.*(u4 + u5*I1).^2 - 2*u1*I2.^2.*(1-u2*I1-u3*I1.^2)+(1-u2*I1-u3*I1.^2).^2);
            end
            
            
        end
    end
end


% Generate 'Check Domain' for the four parts
[CheckDomain_fs] = f_gen_CheckDomain(obj,beamID,IdxEstimation_fs,jT,jGl,nomLoads,pos,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,'FS');
[CheckDomain_rs] = f_gen_CheckDomain(obj,beamID,IdxEstimation_rs,jT,jGl,nomLoads,pos,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,'RS');
[CheckDomain_us] = f_gen_CheckDomain(obj,beamID,IdxEstimation_us,jT,jGl,nomLoads,pos,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,'US');
[CheckDomain_ls] = f_gen_CheckDomain(obj,beamID,IdxEstimation_ls,jT,jGl,nomLoads,pos,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,'LS');

[val_loop,~,ic] = unique([CheckDomain_fs; CheckDomain_rs; CheckDomain_us; CheckDomain_ls],'rows','stable');

beam_elem_loop = cell(size(val_loop,1),1);

% compute exact values of fail. idx at check domain points
for iEl=1:nEl
    for i=1:5
        beam_elem_loop{ic((iEl-1)*5+i)}(end+1) = iEl;
        beam_elem_loop{ic(5*nEl+(iEl-1)*5+i)}(end+1) = iEl;
        beam_elem_loop{ic(10*nEl+(iEl-1)*5+i)}(end+1) = iEl;
        beam_elem_loop{ic(15*nEl+(iEl-1)*5+i)}(end+1) = iEl;
    end
end

obj.str.beam(beamID).calcul_failIdx = 1;

for i=1:size(val_loop,1)
    beam_elem_loop{i} = unique(beam_elem_loop{i});
  
    obj.str.nodal_deflections=disp(val_loop(i,1),:,val_loop(i,2))';
    obj.str.beam(beamID).beam_element_failIdx = beam_elem_loop{i};
    obj.str=obj.str.f_postprocess();
    obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();

    for iEl=1:length(beam_elem_loop{i})
        allFailIdx_fs(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2))=obj.str.beam(beamID).beamelement(beam_elem_loop{i}(iEl)).crosssection.laminate_fs.fail_idx_crit;
        allFailIdx_rs(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2))=obj.str.beam(beamID).beamelement(beam_elem_loop{i}(iEl)).crosssection.laminate_rs.fail_idx_crit;
        allFailIdx_us(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2))=obj.str.beam(beamID).beamelement(beam_elem_loop{i}(iEl)).crosssection.laminate_sk_up.fail_idx_crit;
        allFailIdx_ls(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2))=obj.str.beam(beamID).beamelement(beam_elem_loop{i}(iEl)).crosssection.laminate_sk_lo.fail_idx_crit;
    end
    
end

val_loop_fs = zeros(25*nEl,1);
val_loop_rs = zeros(25*nEl,1);
val_loop_us = zeros(25*nEl,1);
val_loop_ls = zeros(25*nEl,1);

% computation of Interp Domain
[allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,InterpDomain_fs] = f_gen_InterpDomain(obj,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,CheckDomain_fs,'FS',beamID);
[allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,InterpDomain_rs] = f_gen_InterpDomain(obj,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,CheckDomain_rs,'RS',beamID);
[allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,InterpDomain_us] = f_gen_InterpDomain(obj,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,CheckDomain_us,'US',beamID);
[allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,InterpDomain_ls] = f_gen_InterpDomain(obj,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,CheckDomain_ls,'LS',beamID);

% cross check all 4 interp domain so one couple gust length/time step is
% computed only once
for iEl = 1:nEl
    
    for iLoop = 1:5
        val_loop_fs((iEl-1)*25+((iLoop-1)*5+1:iLoop*5),2) = ones(5,1)*InterpDomain_fs{iEl}(iLoop,2);
        val_loop_fs((iEl-1)*25+((iLoop-1)*5+1:iLoop*5),1) = InterpDomain_fs{iEl}(:,1);
        
        val_loop_rs((iEl-1)*25+((iLoop-1)*5+1:iLoop*5),2) = ones(5,1)*InterpDomain_rs{iEl}(iLoop,2);
        val_loop_rs((iEl-1)*25+((iLoop-1)*5+1:iLoop*5),1) = InterpDomain_rs{iEl}(:,1);
        
        val_loop_us((iEl-1)*25+((iLoop-1)*5+1:iLoop*5),2) = ones(5,1)*InterpDomain_us{iEl}(iLoop,2);
        val_loop_us((iEl-1)*25+((iLoop-1)*5+1:iLoop*5),1) = InterpDomain_us{iEl}(:,1);
        
        val_loop_ls((iEl-1)*25+((iLoop-1)*5+1:iLoop*5),2) = ones(5,1)*InterpDomain_ls{iEl}(iLoop,2);
        val_loop_ls((iEl-1)*25+((iLoop-1)*5+1:iLoop*5),1) = InterpDomain_ls{iEl}(:,1);
    end

end

[val_loop,~,ic] = unique([val_loop_fs; val_loop_rs; val_loop_us; val_loop_ls], 'rows');

beam_elem_loop = cell(size(val_loop,1),1);

for iEl=1:nEl
    for i=1:25
        beam_elem_loop{ic((iEl-1)*25+i)}(end+1) = iEl;
        beam_elem_loop{ic(25*nEl+(iEl-1)*25+i)}(end+1) = iEl;
        beam_elem_loop{ic(50*nEl+(iEl-1)*25+i)}(end+1) = iEl;
        beam_elem_loop{ic(75*nEl+(iEl-1)*25+i)}(end+1) = iEl;
    end
end

for i=1:size(val_loop,1)
    
    beam_elem_loop{i} = unique(beam_elem_loop{i});
    
    for iEl=1:length(beam_elem_loop{i})
        if allFailIdx_fs(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2)) == 0
            obj.str.nodal_deflections=disp(val_loop(i,1),:,val_loop(i,2))';
            obj.str.beam(beamID).beam_element_failIdx = beam_elem_loop{i}(iEl);
            obj.str=obj.str.f_postprocess();
            obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
        
            allFailIdx_fs(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2))=obj.str.beam(beamID).beamelement(beam_elem_loop{i}(iEl)).crosssection.laminate_fs.fail_idx_crit;
            allFailIdx_rs(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2))=obj.str.beam(beamID).beamelement(beam_elem_loop{i}(iEl)).crosssection.laminate_rs.fail_idx_crit;
            allFailIdx_us(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2))=obj.str.beam(beamID).beamelement(beam_elem_loop{i}(iEl)).crosssection.laminate_sk_up.fail_idx_crit;
            allFailIdx_ls(val_loop(i,1),beam_elem_loop{i}(iEl),val_loop(i,2))=obj.str.beam(beamID).beamelement(beam_elem_loop{i}(iEl)).crosssection.laminate_sk_lo.fail_idx_crit;
        end
    end
end

%% FINISH COMMENTAIRE ICI


[coeffP_T_fs] = f_interp_InterpDomain(obj,allFailIdx_fs,InterpDomain_fs);
[coeffP_T_rs] = f_interp_InterpDomain(obj,allFailIdx_rs,InterpDomain_rs);
[coeffP_T_us] = f_interp_InterpDomain(obj,allFailIdx_us,InterpDomain_us);
[coeffP_T_ls] = f_interp_InterpDomain(obj,allFailIdx_ls,InterpDomain_ls);

[imax_interpTS_fs] = f_GetMaxInterpDomain(obj,disp,beamID,coeffP_T_fs,InterpDomain_fs,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,'FS');
[imax_interpTS_rs] = f_GetMaxInterpDomain(obj,disp,beamID,coeffP_T_rs,InterpDomain_rs,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,'RS');
[imax_interpTS_us] = f_GetMaxInterpDomain(obj,disp,beamID,coeffP_T_us,InterpDomain_us,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,'US');
[imax_interpTS_ls] = f_GetMaxInterpDomain(obj,disp,beamID,coeffP_T_ls,InterpDomain_ls,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,'LS');

if pos
    if test_maxgusts == 1
        [obj.constraints.gustCases.pos.stress.spFr,varargout{1}] = f_InterpMax(imax_interpTS_fs);
        [obj.constraints.gustCases.pos.stress.spRe,varargout{2}] = f_InterpMax(imax_interpTS_rs);
        [obj.constraints.gustCases.pos.stress.skUp,varargout{3}] = f_InterpMax(imax_interpTS_us);
        [obj.constraints.gustCases.pos.stress.skLo,varargout{4}] = f_InterpMax(imax_interpTS_ls);
    else
        [obj.constraints.gustCases.pos.stress.spFr,~] = f_InterpMax(imax_interpTS_fs);
        [obj.constraints.gustCases.pos.stress.spRe,~] = f_InterpMax(imax_interpTS_rs);
        [obj.constraints.gustCases.pos.stress.skUp,~] = f_InterpMax(imax_interpTS_us);
        [obj.constraints.gustCases.pos.stress.skLo,~] = f_InterpMax(imax_interpTS_ls);
    end
else
    if test_maxgusts == 1
        [obj.constraints.gustCases.neg.stress.spFr,varargout{1}] = f_InterpMax(imax_interpTS_fs);
        [obj.constraints.gustCases.neg.stress.spRe,varargout{2}] = f_InterpMax(imax_interpTS_rs);
        [obj.constraints.gustCases.neg.stress.skUp,varargout{3}] = f_InterpMax(imax_interpTS_us);
        [obj.constraints.gustCases.neg.stress.skLo,varargout{4}] = f_InterpMax(imax_interpTS_ls);
    else
        [obj.constraints.gustCases.neg.stress.spFr,~] = f_InterpMax(imax_interpTS_fs);
        [obj.constraints.gustCases.neg.stress.spRe,~] = f_InterpMax(imax_interpTS_rs);
        [obj.constraints.gustCases.neg.stress.skUp,~] = f_InterpMax(imax_interpTS_us);
        [obj.constraints.gustCases.neg.stress.skLo,~] = f_InterpMax(imax_interpTS_ls);
    end
end

end

function [CheckDomain] = f_gen_CheckDomain(obj,beamID,IdxEstimation,jT,jGl,nomLoads,pos,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,element)
% generate 'Check Domain' : 5 points that should surround the
% global maximum

nT=size(obj.aeSSM.gustData.modes,1);
nL=size(obj.aeSSM.gustData.modes,3);
nEl=length(obj.str.beam(beamID).beamelement);

IdxMax_loads = zeros(nEl,2);
IdxMax_strains = zeros(nEl,2);

sample_time = 1:jT:nT;
sample_gl = 1:jGl:nL;


Loc_Max = cell(nEl,1);
obj.str.beam(beamID).calcul_failIdx = 1;

% computation of all local maximums of the failure index approximation
for iEl = 1:nEl
    Loc_Max{iEl} = [];
    mat_strain = squeeze(IdxEstimation(iEl,sample_time,sample_gl));
    [Max_TS,iGl] = max(mat_strain,[],2);
    
    n_sT = length(sample_time);
    n_sGl = length(sample_gl);
    
    if iGl(1) == n_sGl
        if mat_strain(2,iGl(1)) < Max_TS(1) && mat_strain(2,iGl(1)-1) < Max_TS(1)
            Loc_Max{iEl}(end+1,:) = [1,iGl(1)];
        end
    elseif iGl(1) == 1
        if mat_strain(2,iGl(1)) < Max_TS(1) && mat_strain(2,iGl(1)+1) < Max_TS(1)
            Loc_Max{iEl}(end+1,:) = [1,iGl(1)];
        end
    else
        if mat_strain(2,iGl(1)) < Max_TS(1) && mat_strain(2,iGl(1)+1) < Max_TS(1) && mat_strain(2,iGl(1)-1) < Max_TS(1)
            Loc_Max{iEl}(end+1,:) = [1,iGl(1)];
        end
    end
    
    for iT = 2:n_sT-1
        if iGl(iT) == n_sGl
            if mat_strain(iT+1,iGl(iT)-1) < Max_TS(iT) && mat_strain(iT+1,iGl(iT)) < Max_TS(iT) && mat_strain(iT-1,iGl(iT)-1) < Max_TS(iT) && mat_strain(iT-1,iGl(iT)) < Max_TS(iT)
                Loc_Max{iEl}(end+1,:) = [iT,iGl(iT)];
            end
        elseif iGl(iT) == 1
            if mat_strain(iT+1,iGl(iT)+1) < Max_TS(iT) && mat_strain(iT+1,iGl(iT)) < Max_TS(iT) && mat_strain(iT-1,iGl(iT)+1) < Max_TS(iT) && mat_strain(iT-1,iGl(iT)) < Max_TS(iT)
                Loc_Max{iEl}(end+1,:) = [iT,iGl(iT)];
            end
        else
            if mat_strain(iT+1,iGl(iT)+1) < Max_TS(iT) && mat_strain(iT+1,iGl(iT)) < Max_TS(iT) && mat_strain(iT-1,iGl(iT)+1) < Max_TS(iT) && mat_strain(iT-1,iGl(iT)) < Max_TS(iT)...
                    && mat_strain(iT+1,iGl(iT)-1) < Max_TS(iT) && mat_strain(iT-1,iGl(iT)-1) < Max_TS(iT)
                Loc_Max{iEl}(end+1,:) = [iT,iGl(iT)];
            end
        end
    end
    
    if iGl(n_sT) == n_sGl
        if mat_strain(n_sT-1,iGl(n_sT)) < Max_TS(n_sT) && mat_strain(n_sT-1,iGl(n_sT)-1) < Max_TS(n_sT)
            Loc_Max{iEl}(end+1,:) = [n_sT,iGl(n_sT)];
        end
    elseif iGl(n_sT) == 1
        if mat_strain(n_sT-1,iGl(n_sT)) < Max_TS(n_sT) && mat_strain(n_sT-1,iGl(n_sT)+1) < Max_TS(n_sT)
            Loc_Max{iEl}(end+1,:) = [n_sT,iGl(n_sT)];
        end
    else
        if mat_strain(n_sT-1,iGl(n_sT)) < Max_TS(n_sT) && mat_strain(n_sT-1,iGl(n_sT)+1) < Max_TS(n_sT) && mat_strain(n_sT-1,iGl(n_sT)-1) < Max_TS(n_sT)
            Loc_Max{iEl}(end+1,:) = [n_sT,iGl(n_sT)];
        end
    end
    
    
    Loc_Max{iEl}(:,1) = (Loc_Max{iEl}(:,1)-1)*jT+1;
    Loc_Max{iEl}(:,2) = (Loc_Max{iEl}(:,2)-1)*jGl+1;
    
    Val_Max = zeros(size(Loc_Max{iEl},1),1);
    
    for i=1:size(Loc_Max{iEl},1)
        if allFailIdx_fs(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2)) == 0
            obj.str.nodal_deflections=disp(Loc_Max{iEl}(i,1),:,Loc_Max{iEl}(i,2))';
            obj.str.beam(beamID).beam_element_failIdx = iEl;
            obj.str=obj.str.f_postprocess();
            obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
            
            allFailIdx_fs(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
            allFailIdx_rs(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
            allFailIdx_us(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
            allFailIdx_ls(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
        end
        
        switch element
            case 'FS'
                Val_Max(i) = allFailIdx_fs(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2));
            case 'LS'
                Val_Max(i) = allFailIdx_ls(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2));
            case 'RS'
                Val_Max(i) = allFailIdx_rs(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2));
            case 'US'
                Val_Max(i) = allFailIdx_us(Loc_Max{iEl}(i,1),iEl,Loc_Max{iEl}(i,2));
        end
        
    end
    
    [~,I] = max(Val_Max);
    IdxMax_strains(iEl,:) = Loc_Max{iEl}(I,:);
    
end


% computation of the moments on the wing during gusts
if pos
    for iEl = 1:nEl
 
        Mx = obj.aeSSM.gustData.loads(:,(iEl-1)*6+4,:)+nomLoads((iEl-1)*6+4);
        My = obj.aeSSM.gustData.loads(:,(iEl-1)*6+5,:)+nomLoads((iEl-1)*6+5);
        Mz = obj.aeSSM.gustData.loads(:,(iEl-1)*6+6,:)+nomLoads((iEl-1)*6+6);
        col_mat_l = squeeze(sqrt(Mx.^2+My.^2+Mz.^2));
        
        [~,Il] = max(col_mat_l(:));
        
        [IdxMax_loads(iEl,1),IdxMax_loads(iEl,2)] = ind2sub(size(col_mat_l),Il);
        
    end
else
    for iEl = 1:nEl
        
        Mx = obj.aeSSM.gustData.loads(:,(iEl-1)*6+4,:)-nomLoads((iEl-1)*6+4);
        My = obj.aeSSM.gustData.loads(:,(iEl-1)*6+5,:)-nomLoads((iEl-1)*6+5);
        Mz = obj.aeSSM.gustData.loads(:,(iEl-1)*6+6,:)-nomLoads((iEl-1)*6+6);
        col_mat_l = squeeze(sqrt(Mx.^2+My.^2+Mz.^2));
        [~,Il] = max(col_mat_l(:));
        
        [IdxMax_loads(iEl,1),IdxMax_loads(iEl,2)] = ind2sub(size(col_mat_l),Il);
        
    end
end

IdxMax_loads = IdxMax_strains;

% evaluation of the domain of 5 points based on the moments maximum and the
% failure idx estimation maximum

CheckDomain = zeros(5*nEl,2);
jT = 20; %same issue as for the IdxEstimation
jGl = 2;

for iEl=1:nEl
    if abs(IdxMax_strains(iEl,1) - IdxMax_loads(iEl,1)) <= 1
        if max(IdxMax_strains(iEl,1),IdxMax_loads(iEl,1))+jT > nT
            CheckDomain((iEl-1)*5+[1:5],1) = [nT; nT; min(IdxMax_loads(iEl,1),IdxMax_strains(iEl,1))-jT; min(IdxMax_loads(iEl,1),IdxMax_strains(iEl,1))-jT; nT-jT];
        elseif min(IdxMax_loads(iEl,1),IdxMax_strains(iEl,1))-jT < 1
            CheckDomain((iEl-1)*5+[1:5],1) = [1; 1; 1+2*jT; 1+2*jT; 1+jT];
        else
            CheckDomain((iEl-1)*5+[1:5],1) = [min(IdxMax_loads(iEl,1),IdxMax_strains(iEl,1))-jT; min(IdxMax_loads(iEl,1),IdxMax_strains(iEl,1))-jT;...
                max(IdxMax_strains(iEl,1),IdxMax_loads(iEl,1))+jT; max(IdxMax_strains(iEl,1),IdxMax_loads(iEl,1))+jT; IdxMax_strains(iEl,1)];
        end
    else
        CheckDomain((iEl-1)*5+[1:5],1) = [IdxMax_loads(iEl,1);IdxMax_loads(iEl,1);IdxMax_strains(iEl,1);IdxMax_strains(iEl,1);ceil((IdxMax_loads(iEl,1)+IdxMax_strains(iEl,1))/2)];
    end
    
    if abs(IdxMax_strains(iEl,2) - IdxMax_loads(iEl,2)) <= 1
        if max(IdxMax_strains(iEl,2),IdxMax_loads(iEl,2))+jGl > nL
            CheckDomain((iEl-1)*5+[1:5],2) = [nL;nL-2*jGl;nL;nL-2*jGl;nL-jGl];
        elseif min(IdxMax_loads(iEl,2),IdxMax_strains(iEl,2))-jGl < 1
            CheckDomain((iEl-1)*5+[1:5],2) = [1; 1+2*jGl; 1;1+2*jGl; 1+jGl];
        else
            CheckDomain((iEl-1)*5+[1:5],2) = [min(IdxMax_loads(iEl,2),IdxMax_strains(iEl,2))-jGl;max(IdxMax_strains(iEl,2),IdxMax_loads(iEl,2))+jGl;...
                min(IdxMax_loads(iEl,2),IdxMax_strains(iEl,2))-jGl;max(IdxMax_strains(iEl,2),IdxMax_loads(iEl,2))+jGl;IdxMax_strains(iEl,2)];
        end
    else
        CheckDomain((iEl-1)*5+[1:5],2) = [IdxMax_loads(iEl,2);IdxMax_strains(iEl,2);IdxMax_loads(iEl,2);IdxMax_strains(iEl,2);ceil((IdxMax_loads(iEl,2)+IdxMax_strains(iEl,2))/2)];
    end
end

end

function [allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,InterpDomain] = f_gen_InterpDomain(obj,disp,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,CheckDomain,element,beamID)
%%%% Generate the point for the interpolation given the initial 5 points
%%%% domain computed from loads and strains failure idx estimation

switch element
    case 'FS'
        allFailIdx = allFailIdx_fs;
    case 'RS'
        allFailIdx = allFailIdx_rs;
    case 'LS'
        allFailIdx = allFailIdx_ls;
    case 'US'
        allFailIdx = allFailIdx_us;
end


nT = length(obj.aeSSM.gustData.timeVec);
nL = length(obj.aeSSM.gustData.lengthVec);
nEl = length(squeeze(allFailIdx(1,:,1)));

InterpDomain = cell(nEl,1);

for iEl = 1:nEl
    
    %% Look for the max among the 5 values of CheckDomain
    tst = zeros(5,1);
    
    for i=1:5
        tst(i) = allFailIdx(CheckDomain((iEl-1)*5+i,1),iEl,CheckDomain((iEl-1)*5+i,2));
    end
    [~,I] = max(tst);
    
    InterpDomain{iEl} = zeros(5,2);
    
    %% Get the GL and TS indices of the max and size of CheckDomain
    
    midTS = CheckDomain((iEl-1)*5+I,1);
    gapTS = abs(CheckDomain((iEl-1)*5+5,1)- max(CheckDomain((iEl-1)*5+(1:5),1)));
    
    if gapTS == 1
        gapTS = 2;
    elseif gapTS > 50
        gapTS = ceil(gapTS/4);
    end
    
    midGL = CheckDomain((iEl-1)*5+I,2);
    gapGL = abs(CheckDomain((iEl-1)*5+5,2)- max(CheckDomain((iEl-1)*5+(1:5),2)));
    
    if gapGL == 1
        gapGL = 2;
    elseif gapGL > 15
        gapGL = ceil(gapGL/2);
    end
    
    %% Compute interpolation domain of gust lengths
    
    if midGL == 1 %Check if we are at the border of the domain
        if gapGL > 3
            InterpDomain{iEl}(:,2) = [midGL; midGL + ceil(gapGL/4); midGL + ceil(gapGL/2) ; midGL + ceil(3*gapGL/4); midGL + gapGL];
        else
            InterpDomain{iEl}(:,2) = [midGL; midGL + ceil(gapGL/2); midGL + gapGL ; midGL + ceil(3*gapGL/2); midGL + 2*gapGL];
        end
    elseif midGL == nL
        if gapGL > 3
            InterpDomain{iEl}(:,2) = [midGL - gapGL; midGL - ceil(3*gapGL/4); midGL - ceil(gapGL/2) ; midGL - ceil(gapGL/4); midGL];
        else
            InterpDomain{iEl}(:,2) = [midGL - 2*gapGL; midGL - ceil(3*gapGL/2); midGL - gapGL ; midGL - ceil(gapGL/2); midGL];
        end
        
    else
        
        % Compute a value above and below the reference value to see toward
        % which GL the maximum could be
        val_loop_GL = [midGL - ceil(gapGL/2), midGL + ceil(gapGL/2)];
        
        if val_loop_GL(1) < 1
            val_loop_GL(1) = 1;
        elseif val_loop_GL(2) > nL
            val_loop_GL(2) = nL;
        end
        
        for i=1:length(val_loop_GL)
            if allFailIdx(midTS,iEl,val_loop_GL(i)) == 0

                obj.str.nodal_deflections=disp(midTS,:,val_loop_GL(i))';
                obj.str.beam(beamID).beam_element_failIdx = iEl;
                obj.str=obj.str.f_postprocess();
                obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                
                    allFailIdx_fs(midTS,iEl,val_loop_GL(i))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                    allFailIdx_ls(midTS,iEl,val_loop_GL(i))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                    allFailIdx_rs(midTS,iEl,val_loop_GL(i))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                    allFailIdx_us(midTS,iEl,val_loop_GL(i))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                    
                    switch element
                        case 'FS'
                            allFailIdx(midTS,iEl,val_loop_GL(i)) = allFailIdx_fs(midTS,iEl,val_loop_GL(i));
                        case 'LS'
                            allFailIdx(midTS,iEl,val_loop_GL(i)) = allFailIdx_ls(midTS,iEl,val_loop_GL(i));
                        case 'RS'
                            allFailIdx(midTS,iEl,val_loop_GL(i)) = allFailIdx_rs(midTS,iEl,val_loop_GL(i));
                        case 'US'
                            allFailIdx(midTS,iEl,val_loop_GL(i)) = allFailIdx_us(midTS,iEl,val_loop_GL(i));
                    end
            end
        end
        
        
        
        borne_infGL = val_loop_GL(1);
        borne_supGL = val_loop_GL(2);
        
        [~,imaxGL] = max(squeeze(allFailIdx(midTS,iEl,[borne_infGL midGL borne_supGL])));
        
        if imaxGL == 1 %if the maximum is at the lower indice we need to shift the domain toward smaller value of GL
            if borne_infGL - ceil(gapGL/2) < 1 %always checking that we are staying in the domain
                if gapGL == 1
                    InterpDomain{iEl}(:,2) = [1; midGL; midGL+1; midGL + 2; midGL + 3];
                else
                    InterpDomain{iEl}(:,2) = [1; midGL; midGL+ceil(gapGL/2); midGL + gapGL; midGL + ceil(3*gapGL/2)];
                end
            else
                while borne_infGL - ceil(gapGL/2) > 0 && allFailIdx(midTS,iEl,borne_infGL) > allFailIdx(midTS,iEl,midGL)
                    % check if we are still in the domain before computing
                    % the next value and we have already passed the maximum
                    % or not
                    
                    borne_infGL = borne_infGL - ceil(gapGL/2);
                    
                    if allFailIdx_fs(midTS,iEl,borne_infGL) == 0
                        obj.str.nodal_deflections=disp(midTS,:,borne_infGL)';
                        obj.str.beam(beamID).beam_element_failIdx = iEl;
                        obj.str=obj.str.f_postprocess();
                        obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                        
                        % calculate all failure indices
   
                            allFailIdx_fs(midTS,iEl,borne_infGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                            allFailIdx_ls(midTS,iEl,borne_infGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                            allFailIdx_rs(midTS,iEl,borne_infGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                            allFailIdx_us(midTS,iEl,borne_infGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                            
                            switch element
                                case 'FS'
                                    allFailIdx(midTS,iEl,borne_infGL) = allFailIdx_fs(midTS,iEl,borne_infGL);
                                case 'LS'
                                    allFailIdx(midTS,iEl,borne_infGL) = allFailIdx_ls(midTS,iEl,borne_infGL);
                                case 'RS'
                                    allFailIdx(midTS,iEl,borne_infGL) = allFailIdx_rs(midTS,iEl,borne_infGL);
                                case 'US'
                                    allFailIdx(midTS,iEl,borne_infGL) = allFailIdx_us(midTS,iEl,borne_infGL);
                            end
                    end
                    
                    
                end
                
                if gapGL > 2
                    if borne_infGL ~= midGL - ceil(gapGL/2)
                        InterpDomain{iEl}(:,2) = [borne_infGL; midGL - ceil(gapGL/2); midGL-ceil(gapGL/4); midGL; midGL + ceil(gapGL/2)];
                    else
                        InterpDomain{iEl}(:,2) = [borne_infGL-ceil(gapGL/2); borne_infGL; midGL-ceil(gapGL/4); midGL; midGL + ceil(gapGL/2)];
                    end
                else
                    InterpDomain{iEl}(:,2) = [borne_infGL-ceil(3*gapGL/2); borne_infGL-ceil(gapGL/2); borne_infGL; midGL; midGL + ceil(gapGL/2)];
                end
                
                
            end
            
        elseif imaxGL == 3 %if the maximum is at the bigger indice we need to shift the domain toward bigger value of GL
            
            if  borne_supGL + ceil(gapGL/2) > nL
                if gapGL == 1
                    InterpDomain{iEl}(:,2) = [midGL-3; midGL - 2; midGL-1; midGL;  nL];
                else
                    InterpDomain{iEl}(:,2) = [midGL - ceil(3*gapGL/2); midGL - gapGL; midGL-ceil(gapGL/2); midGL;  nL];
                end
            else
                
                while borne_supGL + ceil(gapGL/2) < nL+1 && allFailIdx(midTS,iEl,borne_supGL) > allFailIdx(midTS,iEl,midGL)
                    
                    borne_supGL = borne_supGL + ceil(gapGL/2);
                    
                    
                    if allFailIdx_fs(midTS,iEl,borne_supGL) == 0
                        obj.str.nodal_deflections=disp(midTS,:,borne_supGL)';
                        obj.str.beam(beamID).beam_element_failIdx = iEl;
                        obj.str=obj.str.f_postprocess();
                        obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                        
                        % calculate all failure indices
                            allFailIdx_fs(midTS,iEl,borne_supGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                            allFailIdx_ls(midTS,iEl,borne_supGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                            allFailIdx_rs(midTS,iEl,borne_supGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                            allFailIdx_us(midTS,iEl,borne_supGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                            
                            switch element
                                case 'FS'
                                    allFailIdx(midTS,iEl,borne_supGL) = allFailIdx_fs(midTS,iEl,borne_supGL);
                                case 'LS'
                                    allFailIdx(midTS,iEl,borne_supGL) = allFailIdx_ls(midTS,iEl,borne_supGL);
                                case 'RS'
                                    allFailIdx(midTS,iEl,borne_supGL) = allFailIdx_rs(midTS,iEl,borne_supGL);
                                case 'US'
                                    allFailIdx(midTS,iEl,borne_supGL) = allFailIdx_us(midTS,iEl,borne_supGL);
                            end
                    end
                    
                    
                    
                end
                
                if gapGL > 2
                    if borne_supGL ~= midGL + ceil(gapGL/2)
                        InterpDomain{iEl}(:,2) = [midGL-ceil(gapGL/2); midGL; midGL + ceil(gapGL/4); midGL + ceil(gapGL/2); borne_supGL];
                    else
                        InterpDomain{iEl}(:,2) = [midGL-ceil(gapGL/2); midGL; midGL + ceil(gapGL/4); borne_supGL; borne_supGL+ceil(gapGL/2)];
                    end
                else
                    InterpDomain{iEl}(:,2) = [midGL - ceil(gapGL/2); midGL; borne_supGL; borne_supGL+ceil(gapGL/2); borne_supGL+ceil(3*gapGL/2)];
                end
                
            end
            
            
            
            
        else %if the maximum was already at the middle value we can directly get a good interpolation domain
            if gapGL > 3
                InterpDomain{iEl}(:,2) = [midGL - ceil(gapGL/2); midGL - ceil(gapGL/4); midGL; midGL + ceil(gapGL/4); midGL + ceil(gapGL/2)];
            else
                InterpDomain{iEl}(:,2) = [midGL - gapGL; midGL - ceil(gapGL/2); midGL; midGL + ceil(gapGL/2); midGL + gapGL];
                
            end
        end
    end
    

    %% Compute the interpolation domain for time steps
    
    if midTS == 1 %Check if we are at the border of the domain
        if gapTS > 3
            InterpDomain{iEl}(:,1) = [midTS; midTS + ceil(gapTS/4); midTS + ceil(gapTS/2) ; midTS + ceil(3*gapTS/4); midTS + gapTS];
        else
            InterpDomain{iEl}(:,1) = [midTS; midTS + ceil(gapTS/2); midTS + gapTS ; midTS + ceil(3*gapTS/2); midTS + 2*gapTS];
        end
    elseif midTS == nT
        if gapTS > 3
            InterpDomain{iEl}(:,1) = [midTS - gapTS; midTS - ceil(3*gapTS/4); midTS - ceil(gapTS/2) ; midTS - ceil(gapTS/4); midTS];
        else
            InterpDomain{iEl}(:,1) = [midTS - 2*gapTS; midTS - ceil(3*gapTS/2); midTS - gapTS ; midTS - ceil(gapTS/2); midTS];
        end
        
    else
        
        % Compute a value above and below the reference value to see toward
        % which GL the maximum could be
        val_loop_TS = [midTS - ceil(gapTS/2), midTS + ceil(gapTS/2)];
        
        
        if val_loop_TS(1) < 1
            val_loop_TS(1) = 1;
        elseif val_loop_TS(2) > nT
            val_loop_TS(2) = nT;
        end
        
        
        for i=1:length(val_loop_TS)
            if allFailIdx_fs(val_loop_TS(i),iEl,midGL) == 0
                obj.str.nodal_deflections=disp(val_loop_TS(i),:,midGL)';
                obj.str.beam(beamID).beam_element_failIdx = iEl;
                obj.str=obj.str.f_postprocess();
                obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                
                % calculate all failure indices

                    allFailIdx_fs(val_loop_TS(i),iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                    allFailIdx_ls(val_loop_TS(i),iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                    allFailIdx_rs(val_loop_TS(i),iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                    allFailIdx_us(val_loop_TS(i),iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                    
                    switch element
                        case 'FS'
                            allFailIdx(val_loop_TS(i),iEl,midGL) = allFailIdx_fs(val_loop_TS(i),iEl,midGL);
                        case 'LS'
                            allFailIdx(val_loop_TS(i),iEl,midGL) = allFailIdx_ls(val_loop_TS(i),iEl,midGL);
                        case 'RS'
                            allFailIdx(val_loop_TS(i),iEl,midGL) = allFailIdx_rs(val_loop_TS(i),iEl,midGL);
                        case 'US'
                            allFailIdx(val_loop_TS(i),iEl,midGL) = allFailIdx_us(val_loop_TS(i),iEl,midGL);
                    end
            end
        end
        
        
        borne_infTS = val_loop_TS(1);
        borne_supTS = val_loop_TS(2);
        
        [~,imaxTS] = max(squeeze(allFailIdx([borne_infTS midTS borne_supTS],iEl,midGL)));
        
        if imaxTS == 1 %if the maximum is at the lower indice we need to shift the domain toward smaller value of GL
            if borne_infTS - ceil(gapTS/2) < 1 %always checking that we are staying in the domain
                              
                if gapTS == 1
                    InterpDomain{iEl}(:,1) = [1; midTS; midTS+1; midTS + 2; midTS + 3];
                else
                    InterpDomain{iEl}(:,1) = [1; midTS; midTS+ceil(gapTS/2); midTS + gapTS; midTS + ceil(3*gapTS/2)];
                end
                
            else
                while borne_infTS - ceil(gapTS/2) > 0 && allFailIdx(borne_infTS,iEl,midGL) > allFailIdx(midTS,iEl,midGL)
                    % check if we are still in the domain before computing
                    % the next value and we have already passed the maximum
                    % or not
                    
                    borne_infTS = borne_infTS - ceil(gapTS/2);
                    
                    
                    
                    if allFailIdx_fs(borne_infTS,iEl,midGL) == 0
                        obj.str.nodal_deflections=disp(borne_infTS,:,midGL)';
                        obj.str.beam(beamID).beam_element_failIdx = iEl;
                        obj.str=obj.str.f_postprocess();
                        obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                        
                        % calculate all failure indices
                            allFailIdx_fs(borne_infTS,iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                            allFailIdx_ls(borne_infTS,iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                            allFailIdx_rs(borne_infTS,iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                            allFailIdx_us(borne_infTS,iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                            
                            switch element
                                case 'FS'
                                    allFailIdx(borne_infTS,iEl,midGL) = allFailIdx_fs(borne_infTS,iEl,midGL);
                                case 'LS'
                                    allFailIdx(borne_infTS,iEl,midGL) = allFailIdx_ls(borne_infTS,iEl,midGL);
                                case 'RS'
                                    allFailIdx(borne_infTS,iEl,midGL) = allFailIdx_rs(borne_infTS,iEl,midGL);
                                case 'US'
                                    allFailIdx(borne_infTS,iEl,midGL) = allFailIdx_us(borne_infTS,iEl,midGL);
                            end
                    end
                    
                    
                    
                end
                
                if gapTS > 3
                    InterpDomain{iEl}(:,1) = [borne_infTS-ceil(gapTS/2); borne_infTS; midTS-ceil(gapTS/4); midTS; midTS + ceil(gapTS/2)];
                else
                    InterpDomain{iEl}(:,1) = [borne_infTS-ceil(3*gapTS/2); borne_infTS-ceil(gapTS/2); borne_infTS; midTS; midTS + ceil(gapTS/2)];
                end
                
            end
            
        elseif imaxTS == 3 %if the maximum is at the bigger indice we need to shift the domain toward bigger value of GL
            
            if borne_supTS + ceil(gapTS/2) > nT+1
               
                if gapTS == 1
                    InterpDomain{iEl}(:,1) = [midTS-3; midTS - 2; midTS-1; midTS;  nT];
                else
                    InterpDomain{iEl}(:,1) = [midTS - ceil(3*gapTS/2); midTS - gapTS; midTS-ceil(gapTS/2); midTS;  nT];
                end
                
            else
                
                while borne_supTS + ceil(gapTS/2) < nT+1 && allFailIdx(borne_supTS,iEl,midGL) > allFailIdx(midTS,iEl,midGL)
                    
                    borne_supTS = borne_supTS + ceil(gapTS/2);
                    
                    
                    if allFailIdx_fs(borne_supTS,iEl,midGL) == 0
                        obj.str.nodal_deflections=disp(borne_supTS,:,midGL)';
                        obj.str.beam(beamID).beam_element_failIdx = iEl;
                        obj.str=obj.str.f_postprocess();
                        obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                        
                        % calculate all failure indices
                        allFailIdx_fs(borne_supTS,iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                        allFailIdx_ls(borne_supTS,iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                        allFailIdx_rs(borne_supTS,iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                        allFailIdx_us(borne_supTS,iEl,midGL)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                        
                        switch element
                            case 'FS'
                                allFailIdx(borne_supTS,iEl,midGL) = allFailIdx_fs(borne_supTS,iEl,midGL);
                            case 'LS'
                                allFailIdx(borne_supTS,iEl,midGL) = allFailIdx_ls(borne_supTS,iEl,midGL);
                            case 'RS'
                                allFailIdx(borne_supTS,iEl,midGL) = allFailIdx_rs(borne_supTS,iEl,midGL);
                            case 'US'
                                allFailIdx(borne_supTS,iEl,midGL) = allFailIdx_us(borne_supTS,iEl,midGL);
                        end
                    end
                    
                    
                    
                end
                if gapTS > 3
                    InterpDomain{iEl}(:,1) = [midTS-ceil(gapTS/2); midTS; midTS + ceil(gapTS/4); borne_supTS; borne_supTS+ceil(gapTS/2)];
                else
                    InterpDomain{iEl}(:,1) = [midTS - ceil(gapTS/2); midTS; borne_supTS; borne_supTS+ceil(gapTS/2); borne_supTS+ceil(3*gapTS/2)];
                end
            end
            
            
            
            
        else %if the maximum was already at the middle value we can directly get a good interpolation domain
            if gapTS > 3
                InterpDomain{iEl}(:,1) = [midTS - ceil(gapTS/2); midTS - ceil(gapTS/4); midTS; midTS + ceil(gapTS/4); midTS + ceil(gapTS/2)];
            else
                InterpDomain{iEl}(:,1) = [midTS - gapTS; midTS - ceil(gapTS/2); midTS; midTS + ceil(gapTS/2); midTS + gapTS];
                
            end
        end
    end
    
    %% Final check to see we are in the domain
    
    if min(InterpDomain{iEl}(:,1)) < 1
        InterpDomain{iEl}(:,1) = InterpDomain{iEl}(:,1) - min(InterpDomain{iEl}(:,1)) + 1;
    elseif max(InterpDomain{iEl}(:,1)) > nT
        InterpDomain{iEl}(:,1) = InterpDomain{iEl}(:,1) - (max(InterpDomain{iEl}(:,1))-nT);
    end
    
    if  min(InterpDomain{iEl}(:,2)) < 1
        InterpDomain{iEl}(:,2) = InterpDomain{iEl}(:,2) - min(InterpDomain{iEl}(:,2)) + 1;
    elseif max(InterpDomain{iEl}(:,2)) > nL
        InterpDomain{iEl}(:,2) = InterpDomain{iEl}(:,2) - (max(InterpDomain{iEl}(:,2))-nL);
    end
    
end

end

function [imax_interpTS] = f_GetMaxInterpDomain(obj,disp,beamID,tot_polT,InterpDomain,allFailIdx_fs,allFailIdx_rs,allFailIdx_us,allFailIdx_ls,element)
        
        switch element
            case 'FS'
                allFailIdx = allFailIdx_fs;
            case 'RS'
                allFailIdx = allFailIdx_rs;
            case 'LS'
                allFailIdx = allFailIdx_ls;
            case 'US'
                allFailIdx = allFailIdx_us;
        end
        
        nEl = length(squeeze(allFailIdx(1,:,1)));
        nT = length(squeeze(allFailIdx(:,1,1)));
        
        imax_interpTS = cell(nEl,1);
        tot_polT_new = tot_polT;
        
        flag_ext_domain = cell(nEl,1);
        
        gapT = zeros(nEl,1);
        imaxT_plot= cell(nEl,1);
        
        for iEl = 1:nEl
            
            flag_ext_domain{iEl} = zeros(length(InterpDomain{iEl}(:,2)),1);
            imaxT_plot{iEl} = zeros(length(InterpDomain{iEl}(:,2)),1);
            imax_interpTS{iEl} = zeros(size(InterpDomain{iEl}(:,2),2),3);
            
            for iGL = 1:length(InterpDomain{iEl}(:,2))
                X = [obj.aeSSM.gustData.timeVec(min(InterpDomain{iEl}(:,1))):0.001:obj.aeSSM.gustData.timeVec(max(InterpDomain{iEl}(:,1)))];
                Z = polyval(fliplr(tot_polT{iEl}(:,iGL)'),X);
                
                
                [maxZ,imaxZ] = max(Z);
                
                tst = abs(obj.aeSSM.gustData.timeVec - X(imaxZ));
                [~, imaxT] = min(tst);
                
                gap = ceil((max(InterpDomain{iEl}(:,1)) - min(InterpDomain{iEl}(:,1)))/2);
                gapT(iEl) = gap;
                
                if imaxZ == length(Z) && max(InterpDomain{iEl}(:,1)) < nT+1
                    
                    if imaxT+gap>nT
                        valT = [imaxT-1 nT];
                    else
                        valT = [imaxT imaxT+gap];
                    end
                    %Calcule combien ca fait en imaxT exactement imaxT + gap/2 ou 4
                    
                    for i=1:length(valT)
                        if allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2)) == 0
                            
                            obj.str.nodal_deflections=disp(valT(i),:,InterpDomain{iEl}(iGL,2))';
                            obj.str.beam(beamID).beam_element_failIdx = iEl;
                            obj.str=obj.str.f_postprocess();
                            obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                            
                            % calculate all failure indices
                            
                            allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                            allFailIdx_rs(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                            allFailIdx_us(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                            allFailIdx_ls(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                            
                            switch element
                                case 'FS'
                                    allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'RS'
                                    allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_rs(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'LS'
                                    allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_ls(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'US'
                                    allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_us(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                            end
                        end
                    end
                    
                    
                    % si la val en gap/2 + 4 est tjrs superieure a celle en imaxT
                    % ben on continue
                    % sinon on prend la valeur en imaxT et imaxT+gap/2 ou4 pour
                    % faire l'interp
                    
                    if allFailIdx(valT(1),iEl,InterpDomain{iEl}(iGL,2)) < allFailIdx(valT(2),iEl,InterpDomain{iEl}(iGL,2))
                        while valT(2)+gap < nT+1 && allFailIdx(valT(1),iEl,InterpDomain{iEl}(iGL,2)) < allFailIdx(valT(2),iEl,InterpDomain{iEl}(iGL,2))
                            valT = [valT(2) valT(2)+gap];
                            for i=1:length(valT)
                                if allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2)) == 0
                                    
                                    obj.str.nodal_deflections=disp(valT(i),:,InterpDomain{iEl}(iGL,2))';
                                    obj.str.beam(beamID).beam_element_failIdx = iEl;
                                    obj.str=obj.str.f_postprocess();
                                    obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                                    
                                    % calculate all failure indices
                                    
                                    allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                                    allFailIdx_rs(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                                    allFailIdx_us(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                                    allFailIdx_ls(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                                    
                                    switch element
                                        case 'FS'
                                            allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                        case 'RS'
                                            allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_rs(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                        case 'LS'
                                            allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_ls(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                        case 'US'
                                            allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_us(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                    end
                                end
                            end
                        end
                        
                        newTS = [InterpDomain{iEl}(1:3) valT(1) valT(2)];
                        
                    else
                        newTS = [InterpDomain{iEl}(1:3) valT(1) valT(2)];
                    end
                    
                    for i=1:length(newTS)
                        if allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) == 0
                            
                            obj.str.nodal_deflections=disp(newTS(i),:,InterpDomain{iEl}(iGL,2))';
                            obj.str.beam(beamID).beam_element_failIdx = iEl;
                            obj.str=obj.str.f_postprocess();
                            obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                            
                            % calculate all failure indices
                            
                            allFailIdx_fs(newTS(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                            allFailIdx_rs(newTS(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                            allFailIdx_us(newTS(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                            allFailIdx_ls(newTS(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                            
                            switch element
                                case 'FS'
                                    allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_fs(newTS(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'RS'
                                    allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_rs(newTS(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'LS'
                                    allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_ls(newTS(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'US'
                                    allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_us(newTS(i),iEl,InterpDomain{iEl}(iGL,2));
                            end
                        end
                    end
                    
                    
                    flag_ext_domain{iEl}(iGL) = 2;
                    
                    imaxT_plot{iEl}(iGL) = valT(2);
                    
                    X = obj.aeSSM.gustData.timeVec(newTS);
                    M = fliplr(vander(X));
                    Ypol = allFailIdx(newTS,iEl,InterpDomain{iEl}(iGL,2));
                    newCoeffPTS = M\Ypol;
                    
                    Xnew = [obj.aeSSM.gustData.timeVec(min(InterpDomain{iEl}(:,1))):0.001:obj.aeSSM.gustData.timeVec(valT(2))];
                    Zext = polyval(fliplr(newCoeffPTS'),Xnew);
                    [maxZ, imaxZ] = max(Zext);
                    
                    imax_interpTS{iEl}(iGL,:) = [Xnew(imaxZ),obj.aeSSM.gustData.lengthVec(InterpDomain{iEl}(iGL,2)),maxZ];
                    tot_polT_new{iEl}(:,iGL) = newCoeffPTS;
                    
                elseif imaxZ == 1 && min(InterpDomain{iEl}(:,1)) > 0
                    flag_ext_domain{iEl}(iGL) = 1;
                    
                    if imaxT-gap < 1
                        valT = [1 imaxT+1];
                    else
                        valT = [imaxT-gap imaxT];
                    end
                    
                    %Calcule combien ca fait en imaxT exactement imaxT + gap/2 ou 4
                    
                    for i=1:length(valT)
                        if allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2)) == 0
                            obj.str.nodal_deflections=disp(valT(i),:,InterpDomain{iEl}(iGL,2))';
                            obj.str.beam(beamID).beam_element_failIdx = iEl;
                            obj.str=obj.str.f_postprocess();
                            obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                            
                            % calculate all failure indices
                            allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                            allFailIdx_rs(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                            allFailIdx_us(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                            allFailIdx_ls(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                            
                            switch element
                                case 'FS'
                                    allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'RS'
                                    allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_rs(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'LS'
                                    allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_ls(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'US'
                                    allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_us(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                            end
                            
                        end
                    end
                    
                    % si la val en gap/2 + 4 est tjrs superieure a celle en imaxT
                    % ben on continue
                    % sinon on prend la valeur en imaxT et imaxT+gap/2 ou4 pour
                    % faire l'interp
                    
                    if allFailIdx(valT(1),iEl,InterpDomain{iEl}(iGL,2)) > allFailIdx(valT(2),iEl,InterpDomain{iEl}(iGL,2))
                        while valT(1)-gap > 1 && allFailIdx(valT(1),iEl,InterpDomain{iEl}(iGL,2)) > allFailIdx(valT(2),iEl,InterpDomain{iEl}(iGL,2))
                            valT = [valT(1)-gap valT(1)];
                            
                            for i=1:length(valT)
                                if allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) == 0
                                    obj.str.nodal_deflections=disp(valT(i),:,InterpDomain{iEl}(iGL,2))';
                                    obj.str.beam(beamID).beam_element_failIdx = iEl;
                                    obj.str=obj.str.f_postprocess();
                                    obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                                    
                                    % calculate all failure indices
                                    allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                                    allFailIdx_rs(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                                    allFailIdx_us(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                                    allFailIdx_ls(valT(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                                    
                                    switch element
                                        case 'FS'
                                            allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_fs(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                        case 'RS'
                                            allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_rs(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                        case 'LS'
                                            allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_ls(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                        case 'US'
                                            allFailIdx(valT(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_us(valT(i),iEl,InterpDomain{iEl}(iGL,2));
                                    end
                                end
                            end
                            
                            
                        end
                        
                        newTS = [valT(1) valT(2) InterpDomain{iEl}(3:5)];
                        
                    else
                        newTS = [valT(1) valT(2) InterpDomain{iEl}(3:5)];
                    end
                    
                    for i=1:length(newTS)
                        if allFailIdx_fs(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) == 0
                            obj.str.nodal_deflections=disp(newTS(i),:,InterpDomain{iEl}(iGL,2))';
                            obj.str.beam(beamID).beam_element_failIdx = iEl;
                            obj.str=obj.str.f_postprocess();
                            obj.str.beam(beamID) = obj.str.beam(beamID).f_calc_stresses();
                            
                            % calculate all failure indices
                            
                            allFailIdx_fs(newTS(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                            allFailIdx_rs(newTS(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                            allFailIdx_us(newTS(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                            allFailIdx_ls(newTS(i),iEl,InterpDomain{iEl}(iGL,2))=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                            
                            switch element
                                case 'FS'
                                    allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_fs(newTS(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'RS'
                                    allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_rs(newTS(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'LS'
                                    allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_ls(newTS(i),iEl,InterpDomain{iEl}(iGL,2));
                                case 'US'
                                    allFailIdx(newTS(i),iEl,InterpDomain{iEl}(iGL,2)) = allFailIdx_us(newTS(i),iEl,InterpDomain{iEl}(iGL,2));
                            end
                        end
                    end
                    
                    imaxT_plot{iEl}(iGL) = valT(1);
                    
                    X = obj.aeSSM.gustData.timeVec(newTS);
                    M = fliplr(vander(X));
                    Ypol = allFailIdx(newTS,iEl,InterpDomain{iEl}(iGL,2));
                    newCoeffPTS = M\Ypol;
                    
                    Xnew = [obj.aeSSM.gustData.timeVec(valT(1)):0.001:obj.aeSSM.gustData.timeVec(max(InterpDomain{iEl}))];
                    Zext = polyval(fliplr(newCoeffPTS'),Xnew);
                    [maxZ, imaxZ] = max(Zext);
                    
                    imax_interpTS{iEl}(iGL,:) = [Xnew(imaxZ),obj.aeSSM.gustData.lengthVec(InterpDomain{iEl}(iGL,2)),maxZ];
                    tot_polT_new{iEl}(:,iGL) = newCoeffPTS;
                    
                else
                    imax_interpTS{iEl}(iGL,:) = [X(imaxZ),obj.aeSSM.gustData.lengthVec(InterpDomain{iEl}(iGL,2)),maxZ];
                end
                
            end
            
        end
        
    end

function [coeffP_T] = f_interp_InterpDomain(obj,allFailIdx,InterpDomain)

nEl = length(squeeze(allFailIdx(1,:,1))); 

coeffP_T = cell(nEl,1);

for iEl = 1:nEl
    coeffP_T{iEl} = zeros(5,5);
    
    X = obj.aeSSM.gustData.timeVec(InterpDomain{iEl}(:,1));
    M = fliplr(vander(X));
    
    for iGL = 1:size(InterpDomain{iEl}(:,2),1)
        Ypol = allFailIdx(InterpDomain{iEl}(:,1),iEl,InterpDomain{iEl}(iGL,2));
        coeffP_T{iEl}(:,iGL) = M\Ypol;
    end
       
end
end

function [interp_max_out,interp_max] = f_InterpMax(imax_interpTS)

nEl = length(imax_interpTS);

interp_max = zeros(nEl,3);

a = zeros(nEl,1);
b = zeros(nEl,1);

coeffPolPl = cell(nEl,1);

for iEl = 1:nEl
    Y = imax_interpTS{iEl}(:,2);
    X = [ones(length(Y),1) imax_interpTS{iEl}(:,1)];
    
    if max(X(:,2))-min(X(:,2)) == 0
        Xpl = imax_interpTS{iEl}(:,2);
        M = fliplr(vander(Xpl));
        Ypl = imax_interpTS{iEl}(:,3);
        
        coeffPolPl{iEl} = pinv(M)*Ypl;
        %pinv used here because sometimes the values in M varies from order 10^3 to order 10^9 and more which leads to a badly scaler matrix that Matlab does not like and issue a warning for
        
        Xpol = min(Xpl):0.001:max(Xpl);
        Ypol = polyval(fliplr(coeffPolPl{iEl}'),Xpol);
        
        [Zmax,J] = max(Ypol);
        
        Xmax = imax_interpTS{iEl}(1,1);
        Ymax = Xpol(J);
    else
        B = X\Y;
        a(iEl) = B(2);
        b(iEl) = B(1);

        [~, I] = sort(imax_interpTS{iEl}(:,3),'descend');
        val_5max = imax_interpTS{iEl}(I(1:5),:);
        
        Xpl = sqrt(val_5max(:,1).^2+(val_5max(:,2)-b(iEl)).^2)/10;
        M = fliplr(vander(Xpl));
        Ypl = val_5max(:,3);
        
        coeffPolPl{iEl} = pinv(M)*Ypl;
        %pinv used here because sometimes the values in M varies from order 10^3 to order 10^9 and more which leads to a badly scaler matrix that Matlab does not like and issue a warning for
        
        Xpol = min(Xpl):0.001:max(Xpl);
        Ypol = polyval(fliplr(coeffPolPl{iEl}'),Xpol);
        
        [Zmax,J] = max(Ypol);
        
        Xmax = 10*Xpol(J)/(sqrt(1+a(iEl)^2));
        Ymax = a(iEl)*Xmax+b(iEl);
    end
    
    interp_max(iEl,:) = [Xmax Ymax Zmax];
end

interp_max_out = interp_max(:,3)';

end