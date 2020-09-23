function entropy = calc_entropy( spatial, info_len)
        rho = HILL(spatial);
%         rho = MiPOD(spatial,min(info_len,round(.7*numel(spatial)))./numel(spatial));

        rho0 = rho*0;
        wetConst = 1e13;
        rhoM1 = rho;
        rhoP1 = rho;
        rhoP1(rhoP1 > wetConst) = wetConst; 
        rhoP1(spatial == 255) = wetConst;

        rhoM1(rhoM1 > wetConst) = wetConst;
        rhoM1(spatial == 0) = wetConst;

        [~, ~, pM, p0, pP]  =  EmbeddingSimulator_hu( spatial, rhoM1, rho0, rhoP1, min(info_len,round(.7*numel(spatial))));        
        pM(pM<1e-7)=1e-7;
        pP(pP<1e-7)=1e-7;
        entropy = -pM.*log2(pM) - p0.*log2(p0) - pP.*log2(pP);
end

