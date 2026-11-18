function plot_results(frequencies, Gamma, T)
    freq_MHz = frequencies / 1e6;
    mag_Gamma = abs(Gamma);   mag_T = abs(T);
    dB_Gamma  = 20*log10(mag_Gamma);
    dB_T      = 20*log10(mag_T);

    figure;
    subplot(2,1,1);
    plot_pair(freq_MHz, mag_Gamma, mag_T, @max, ...
              'Magnitude', 'Electric Reflection and Transmission Coefficients (Relative)', ...
              '%.3f', '');

    subplot(2,1,2);
    plot_pair(freq_MHz, dB_Gamma, dB_T, @max, ...
              'Magnitude (dB)', 'Electric Reflection and Transmission Coefficients (dB)', ...
              '%.2f', ' dB');
end

function plot_pair(f, yG, yT, peakFcn, ylab, ttl, fmt, unit)
    [mG, iG] = peakFcn(yG);
    [mT, iT] = peakFcn(yT);

    plot(f, yG, 'b', 'LineWidth', 1.5); hold on;
    plot(f, yT, 'r', 'LineWidth', 1.5);
    plot(f(iG), mG, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    plot(f(iT), mT, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    hold off;

    xlabel('Frequency (MHz)'); ylabel(ylab); title(ttl);
    legend('|\Gamma|', '|T|', ...
           sprintf(['Max |\\Gamma| = ' fmt unit ' @ %.1f MHz'], mG, f(iG)), ...
           sprintf(['Max |T| = '      fmt unit ' @ %.1f MHz'], mT, f(iT)), ...
           'Location', 'best');
    grid on; xlim([min(f), max(f)]);
end