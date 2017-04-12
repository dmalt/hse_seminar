function [A, Ps, Cs, IND] = DICS(C, G, lambda, is_imag)
% -------------------------------------------------------
% DICS beamformer on cross-spectrum. compute beamformer
% topographies and use them to estimate powers and 
% rank source-level connections by their correlation
% with sensor-level cross-spectrum
% -------------------------------------------------------
% FORMAT:
%  [A, Ps, Cs, IND] = DICS(C, G, lambda, is_imag)
% INPUTS:
%   C        - {n_sensors x n_sensors} matrix;
%              sensor-level cross-spectrum
%   G        - {n_sensors x n_sources * 2} matrix;
%              forward operator
%   lambda   - scalar; matrix shrinkage parameter
%              default = 10;
%   is_imag  - boolean flag; if true perform iDICS version
%              of the algorithm
% OUTPUTS:
%   A        - {n_sources * 2 x n_sensors} matrix;
%              inverse operator
%   Ps       - {n_sensors x 1} vector;
%              source-level power estimates
%   Cs       - {n_sources * (n_sources - 1) / 2 x 1} vector
%              vectorized upper diagonal of matrix of
%              source-level connection correlations with
%              sensor-level cross-spectrum
%   IND      - {n_sources * (n_sources - 1) / 2 x 2} matrix;
%              mapping between linear and matrix index notation;
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    import ps.PSIICOS_ScanFast

    if nargin < 4
        is_imag = false;
    end

    if nargin < 3
        lambda = 10;
    end

    Nch = size(C,1);
    C_reg = inv(C + lambda * trace(C) / Nch * eye(Nch));
    Ns = fix(0.5 * size(G, 2)); % assume tangent space dimension of 2

    range = 1:2;
    A = zeros(size(G'));
    % A = G';

    for i = 1:Ns
        L = G(:,range);
        A(range,:) = inv(L' * C_reg * L) * L' * C_reg;
        range = range + 2;
    end

    % ---------- Estimate powers -------- %
    Ps = zeros(Ns,1);
    for i=1:Ns
        range_i = i * 2 - 1 : i * 2;
        ai = A(range_i,:);
        cs = ai * C * ai';
        [~, s, ~] = svd(cs);
        Ps(i) = sqrt(s(1,1));
    end;
    % ----------------------------------- %

    % Rank connections by their correlation with vec(C)
    [Cs, IND] = FastLambda1(A', C(:), is_imag);
end

function [Cs, IND] = FastLambda1(G2dU, Cp, is_imag)
% -------------------------------------------------------------------------
% Find correlations of cross-spectrum with the forward operator
% for each cortical source
% -------------------------------------------------------------------------
% FORMAT:
%   [Cs, IND] = FastLambda1(G2dU, Cp, is_imag) 
% INPUTS:
%   G2dU       - {n_sensors x n_sources} matrix of forward model
%   Cp         - {n_sensors ^ 2 x n_times} or
%                {n_sensors ^ 2 x n_components}
%                matrix of timeseries or left singular vectors of
%                cross-spectrum on sensors
%   is_imag    - boolean flag; if true perform imaginary version of
%                algorithm; default = false
% OUTPUTS:
%   Cs         - {(n_sources ^ 2 - n_sources) / 2} matrix of 
%                correlations between source topographies
%                and forward operator
%   IND        - {(n_sources ^ 2 - n_sources) / 2} matrix of
%                indices to build a mapping between upper
%                triangle and (i,j) matrix indexing 
%                IND(l,:) --> [i,j]
% ___________________________________________________________________________
% Alex Ossadtchii, ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    if nargin < 3
        is_imag = false;
    end

    [Nsns, Nsrc2] = size(G2dU);
    Nsrc = Nsrc2 / 2;

    % if(size(Cp,1)~= Nsns)
    %     disp('Incomptible dimensions G2dU vs Cp');
    %     return;
    % end

    n_comp = size(Cp, 2);
    T = zeros(n_comp, Nsrc * (Nsrc - 1) / 2);
    D = zeros(n_comp, Nsrc * (Nsrc - 1) / 2);
    IND = zeros(Nsrc * (Nsrc - 1) / 2, 2);

    % below is the optimized implementation of this:
    % Look at each pair and estimate subspace correlation
    % between cross-spectrum and topography of this pair
    % tic
    % p = 1;
    % for i=1:Nsrc
    %     range_i = i * 2 - 1 : i * 2;
    %     ai = G2dU(:,range_i)';
    %     for j=i + 1:Nsrc
    %          range_j = j * 2 - 1 : j * 2;
    %          aj = G2dU(:, range_j)';
    %         cs = ai * Cp * aj';
    %         [u s v] = svd(cs);
    %         Cs(p) = max(diag(s)); 
    %          p = p + 1;
    %     end;
    % end;
    % toc
    %

    p = 1;
    for iSrc = 1:Nsrc
        for iComp = 1:n_comp
            Cp_sq = reshape(Cp(:,iComp), Nsns, Nsns);
            % --- Take iSrc-th location topographies ---- %
            ai = G2dU(:, iSrc * 2 - 1 : iSrc * 2)';       
            tmp = ai * Cp_sq;
            cslong = tmp * G2dU;
            if is_imag
                cslong = imag(cslong);
            end
            cs2long = cslong .* conj(cslong);
            cs2longd = cslong(1,:) .* conj(cslong(2,:));
            cs2_11_22 = [sum(reshape(cs2long(1,:), 2, Nsrc), 1);...
                         sum(reshape(cs2long(2,:), 2, Nsrc), 1)];
            cs2_12_21 = sum(reshape(cs2longd, 2, Nsrc), 1);
            Ti = sum(cs2_11_22, 1);
            Di = prod(cs2_11_22, 1) - cs2_12_21 .* conj(cs2_12_21);
            T(iComp, p : p + Nsrc - 1 - iSrc) = Ti(iSrc + 1 : Nsrc);
            D(iComp, p : p + Nsrc - 1 - iSrc) = Di(iSrc + 1 : Nsrc);
        end
        IND(p : p + Nsrc - iSrc - 1, 2) = (iSrc + 1 : Nsrc)';
        IND(p : p + Nsrc - iSrc - 1, 1) = iSrc;
        p = p + Nsrc - iSrc;
    end;
    Cs = 0.5 * T + sqrt(0.25 * T .* T - D); 
    % Cs = sum(Cs,1);  
    Cs = max(Cs,[],1);    
end
