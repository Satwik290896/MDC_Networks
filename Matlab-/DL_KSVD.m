function mD = DL_KSVD(mInitD, mX, cardinality)
      mD            = mInitD;
      numberOfAtoms = size(mD, 2);
      for ii = 1 : 50
          %-- Update the Representations:
          mA = omp(mX, mD, cardinality);
          %-- Update the Dictionary:
          for jj = 1 : 5
              mE = mX - mD * mA;
              for kk = 1 : numberOfAtoms
                  vP        = find(mA(kk,:));
                  mEP       = mE(:, vP) + mD(:,kk) * mA(kk, vP);
                  vA        = mD(:,kk)' * mEP;
                  mA(kk,vP) = vA;
                  mD(:,kk)  = mEP * vA' / (vA * vA');
              end
          end
          mD = bsxfun(@rdivide, mD, sqrt(sum(mD.^2, 1)));
          %% Progress (Debug):
  %         mA     = omp(mD' * mX, mD' * mD, cardinality);
  %         mX_hat = mD * mA;
  %         A = mean(abs(mX(:) - mX_hat(:)))
      end
end