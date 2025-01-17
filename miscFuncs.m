classdef miscFuncs

    properties
    end

    methods (Static)

        function normSignal = minmaxNorm(signal)
            % assumes rows to be trials

            sMax = max(signal,[],2);
            sMin = min(signal, [],2);

            normSignal = ((signal - sMin)./(sMax - sMin));
        end

        function normSignal = minmaxNorm_Minus1ToPlus1(signal)
            % assumes rows to be trials

            sMax = max(signal,[],2);
            sMin = min(signal, [],2);

            normSignal = ((signal - sMin)./(sMax - sMin))*2 -1;
        end


        function animateplot(x,y,step_size, xlab,ylab,linecolor,linewidth,markercolor,saveflag)

            figure("WindowState", "maximized");
            pause(1e-6)
            plot(x,y,'Color','none');
            hold on;

            p = plot(x,y,'Color',linecolor,'LineWidth',linewidth);
            m = scatter(x, y, 'filled','o','MarkerFaceColor',markercolor, 'LineWidth',1);
            xlabel(xlab);
            ylabel(ylab);
            a= gca; a.XAxis.Visible="off"; a.YAxis.Visible="off";
            box off;

            for k = 1:step_size:length(x)

                p.XData = x(1:k);
                p.YData = y(1:k);

                m.XData = x(k);
                m.YData = y(k);

                pause(1e-6)



                % Saving the figure
                if saveflag==1
                    filename = 'animation.gif';
                    frame = getframe(gcf);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    if k == 1
                        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
                            'DelayTime',0.1);
                    else
                        imwrite(imind,cm,filename,'gif','WriteMode','append',...
                            'DelayTime',0.1);
                    end
                end
            end

        end

        function [pvalA, pvalAB] = returnPvals(A,B)
            % Calculates p-values of all pairs of columns in matrix A and B
            % pvalue = ranksum(A(:,i),B(:,i))
            % B is assumed to be baseline values, not testing them pairwise

            ncol = size(A,2);
            pvalA = zeros([ncol,ncol]);
            pvalAB = zeros([1,ncol]);

            for i = 1:ncol
                for j = 1:ncol
                    pvalA(i,j) = signrank(A(:,i), A(:,j));
                end

                if ~isempty(B)
                    pvalAB(i) = signrank(A(:,i),B(:,i));
                else
                    pvalAB = [];
                    continue;
                end
            end
        end

        function [onLoc, offLoc] = findSSbounds(signal, prcThresh, nPts, fs)

            narginchk(1,4)
            if nargin < 2
                prcThresh = 0.95;
                nPts = 50;
            elseif nargin < 3
                nPts = 50;
            end
            signal = abs(signal);
            % figure;
            % plot(signal); hold on;

            if length(signal) < fs+1
                if mod(length(signal),2) == 0
                    frameLength = length(signal)-1;
                else
                    frameLength = length(signal);
                end
            else
                frameLength = fs+1;
            end

            signal = sgolayfilt(signal, 3, frameLength);
            threshold = prcThresh*max(signal);
            % plot(signal);
            % yline(threshold, 'k--');
            risingLocs = zeros([1 length(signal)]);
            fallinsLocs = zeros([1 length(signal)]);

            for isample = (nPts+1):(length(signal)-nPts)

                if mean(signal((isample-nPts):isample)) <  threshold && mean(signal((isample+1):(isample+nPts))) > threshold
                    risingLocs(isample) = 1;
                elseif mean(signal((isample-nPts):isample)) >  threshold && mean(signal((isample+1):(isample+nPts))) < threshold
                    fallinsLocs(isample) = 1;
                else
                    continue;
                end
            end

            onLoc = find(risingLocs, 1, "first" );
            offLoc = find(fallinsLocs, 1, "last" );

            % xline(onLoc, 'g--');
            % xline(offLoc, 'r--');
        end

        function [onLoc, offLoc] = threshCrossIdx(signal, thresh, nPts, orderFlag)


            signal = abs(signal);
            threshold = thresh*max(signal);
            % figure;
            % plot(stim); hold on;
            % yline(threshold, 'k--');
            risingLocs = zeros([1 length(signal)]);
            fallingLocs = zeros([1 length(signal)]);
            signalFlipped = flip(signal);
            for isample = (nPts+1):(length(signal)-nPts)

                if mean(signal((isample-nPts):isample)) <  threshold && mean(signal((isample+1):(isample+nPts))) > threshold
                    risingLocs(isample) = 1;
                end
                % if mean(signal((isample-nPts):isample)) >  threshold && mean(signal((isample+1):(isample+nPts))) < threshold
                if mean(signalFlipped((isample-nPts):isample)) <  threshold && mean(signalFlipped((isample+1):(isample+nPts))) > threshold
                    fallingLocs(isample) = 1;
                end
            end
            fallingLocs = flip(fallingLocs);
            onLoc = find(risingLocs, 1, 'first' );
            length(find(fallingLocs))
            offLoc = find(fallingLocs,1, "first");

            % xline(find(risingLocs, 1 ), 'g--');
            % xline(find(fallinsLocs, 1,"last" ), 'r--');
        end

        function out_range = range_converter(in_range, anchor, zoom)

            in_axis_range = diff(in_range);
            left = (anchor - in_range(1)) / in_axis_range;
            right = (in_range(2) - anchor) / in_axis_range;
            out_axis_range = in_axis_range / zoom;
            out_range = [anchor - left * out_axis_range, anchor + right * out_axis_range];
        end

        function vel = getVelocity(stim, fs)
            [b,a] = butter(3,4/(fs/2), 'low');
            velocity = diff(stim)*fs;
            vel_filtered = filtfilt(b,a,velocity);
            vel = max(vel_filtered);
        end

        function [low, high, med] = getWhisker(data)

        iqr_data = iqr(data);

        low = prctile(data,25) - 1.5*iqr_data;
        high = prctile(data, 75) + 1.5*iqr_data;
        med = median(data);

        end

        function f = getFreqFit(stim, fs, ON_dur)

            stim_smooth = sgolayfilt(stim,2,101);
            [c, midlev] = midcross(stim_smooth, fs, 'Tolerance',2);
            t = linspace(0,ON_dur, length(stim));
            % figure;
            % plot(t, stim_smooth); hold on;
            % plot(c',midlev,'rx')
            
            c = [0 c];
            T_half = diff(c);
            freq = 1./(2*T_half);
            t_f = T_half./2;
            t_f = c(1:end-1)+T_half./2;
            
            % figure; plot(t,stim_smooth); xline(t_f, 'k--');
            
            f = fit(t_f',freq','poly1');
            % figure; plot(f,t_f,freq);

        end

        function [p_wt, f_wt] = getPSD(signal, fs, FreqLims)
            
            [wt, f_wt] = cwt(signal, "amor", fs, VoicesPerOctave=48, FrequencyLimits = FreqLims);

            % figure;
            % cwt(signal, "amor", fs);
            
            p_wt = mean(abs(wt), 2);
        end

    end

    
end



