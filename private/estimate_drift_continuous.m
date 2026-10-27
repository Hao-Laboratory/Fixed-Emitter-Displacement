function [d_est, se_est, keep] = estimate_drift_continuous(dispXY, madTau)
%ESTIMATE_DRIFT_CONTINUOUS Mean [dx dy], component-wise SE, and MAD inlier mask.
% madTau <= 0 disables rejection; a zero MAD scale disables that axis.

    if isempty(dispXY)
        d_est = [nan nan];
        se_est = [nan nan];
        keep = false(0,1);
        return;
    end

    dx = dispXY(:,1);
    dy = dispXY(:,2);
    keep = isfinite(dx) & isfinite(dy);

    if madTau > 0 && any(keep)
        mx = median(dx(keep));
        my = median(dy(keep));
        sx = 1.4826 * median(abs(dx(keep) - mx));
        sy = 1.4826 * median(abs(dy(keep) - my));

        if ~(isfinite(sx) && sx > 0), sx = inf; end
        if ~(isfinite(sy) && sy > 0), sy = inf; end

        keep = keep & abs(dx - mx) <= madTau * sx & abs(dy - my) <= madTau * sy;
    end

    if ~any(keep)
        d_est = [nan nan];
        se_est = [nan nan];
        return;
    end

    dxk = dx(keep);
    dyk = dy(keep);
    K = numel(dxk);

    d_est = [mean(dxk), mean(dyk)];
    se_est = [sqrt(var(dxk,1)/K), sqrt(var(dyk,1)/K)];
end
