classdef miscFuncs

    properties
    end

    methods (Static)

        function normSignal = minmaxNorm(signal)
            % assumes rows to be trials

            sMax = max(signal,[],2);
            sMin = min(signal, [],2);

            normSignal = (signal - sMin)./(sMax - sMin);
        end


        function animateplot(x,y,xlab,ylab,linecolor,linewidth,markercolor,saveflag)


            plot(x,y,'Color','none');
            hold on;

            p = plot(x,y,'Color',linecolor,'LineWidth',linewidth);
            m = scatter(x, y, 'filled',markercolor);
            xlabel(xlab);
            ylabel(ylab);


            for k = 1:20:length(x)

                p.XData = x(1:k);
                p.YData = y(1:k);

                m.XData = x(k);
                m.YData = y(k);

                pause(1e-6)


                % Saving the figure
                if saveflag==1
                    filename = 'animation';
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
                    pvalA(i,j) = ranksum(A(:,i), A(:,j));
                end

                pvalAB(i) = ranksum(A(:,i),B(:,i));
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
    end
end



