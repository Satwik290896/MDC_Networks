function [mA] = DL_KSVD_test(mInitD, mX)
      mD            = mInitD;
      numberOfAtoms = size(mD, 2);
      for ii = 1 : 1
          %-- Update the Representations:
          %mA = omp(mX, mD, cardinality);
          opt = struct('ReguParaType','SLEP','ReguAlpha',0.05)
          mA = SparseCodingwithBasis(mD,mX,opt)
          mA = cell2mat(mA)';
          %-- Update the Dictionary:
          %for jj = 1 : 5
          %    mE = mX - mD * mA;
          %    for kk = 1 : numberOfAtoms
          %        vP        = find(mA(kk,:));
          %        mEP       = mE(:, vP) + mD(:,kk) * mA(kk, vP);
          %        vA        = mD(:,kk)' * mEP;
          %        mA(kk,vP) = vA;
          %        mD(:,kk)  = mEP * vA'%/ (vA * vA');
          %    end
          %end
          %mD = bsxfun(@rdivide, mD, sqrt(sum(mD.^2, 1)));
          %% Progress (Debug):
  %         mA     = omp(mD' * mX, mD' * mD, cardinality);
  %         mX_hat = mD * mA;
  %         A = mean(abs(mX(:) - mX_hat(:)))
      end
end