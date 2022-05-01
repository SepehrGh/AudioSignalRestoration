function [lowThreshold, highThreshold, clippedIndexesLow, clippedIndexesHigh] = extractClippingLevel(signal)
    lowThreshold = min(signal);
    highThreshold = max(signal);
    clippedIndexesLow  = (signal == highThreshold);
    clippedIndexesHigh = (signal == lowThreshold);
end
%     filename_list = ["a08_violin", "a16_clarinet", "a18_bassoon", "a25_harp", "a35_glockenspiel", ...
%                      "a41_celesta", "a42_accordion", "a58_guitar_sarasate", "a60_piano_schubert",  ...
%                      "a66_wind_ensemble_stravinsky"];
%     decibel_list = ["01", "03", "05", "07", "10", "15", "20"];
%     clValues = zeros(10,7);
%     for i=1:length(decibel_list)
%         decibel = decibel_list(i);
%         dirAddr = "C:\Users\Sepehr\Desktop\Audio Restoration\SurveySamples\" + decibel + "dB\";
%         for j=1:length(filename_list)
%             filename = filename_list(j);
%             fileAddr = dirAddr + filename + "_" + decibel + ".wav";
%             [data,Fs] = audioread(fileAddr);
%             clValues(j,i) = max(abs(data));
%         end
%     end
%     resMat = horzcat(transpose(filename_list),clValues);
%     res = vertcat(horzcat("Sample\\dB", decibel_list), resMat)
