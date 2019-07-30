function f = read_all_xvals(prefix, suffix, outputFileName, ...
                                unitsFlag, plotFlag,...
                                dt, tmaxvals, xvals, delay, sigma, offset)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function reads a dumbbell data file and performs a TEXTRA
    % extrapolation on the data from different timestep widths to give
    % results suitable for further MATLAB post-processing
    %
    % direct inputs:
    %   inputFileName - the name of the dumbbell code output data file
    %   outputFileName - the name of the output re-arranged file
    %   unitsFlag - 0 for no changes, 1 for rodlike viscosity data, 2 for
    %       rodlike normal stress difference data, 3 for rodlike time only,
    %       4 for chiTau/chiG data (absolute values)
    %   plotFlag - 0 for no plot, 1 for plot of data
    %   dt - horizontal array of timestep widths for input file
    %   tmaxvals - total simulation runtime array
    %   xvals - independent variable values array
    %   delay - time delay between data output to file
    %   sigma - value of sigma in simulation units
    %   offset - this only exists because of a silly bug I had in previous
    %       code versions which made certain output files slightly
    %       different. It should just be set to 0 now.
    %
    % direct outputs:
    %   f - the same data which is placed in the output file
    % 
    % implicit outputs:
    %   saveFile - file with name outputFileName which contains the TEXTRA
    %       extrapolated data, sorted by time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j=1:length(xvals)
        str = sprintf('%.15f ',xvals(j));
        if floor(xvals(j))==xvals(j)
            str = regexprep(str, '\.[0]+ ', '');
        else
            str = regexprep(str, '[0]+ ', '');
        end
        file = strcat(prefix, str, suffix);
    
        tmax = tmaxvals(j);
        
        rawData = dlmread(file, '');
        size_steps=floor(tmax/delay+1+offset);

        data = nan(size_steps,3,length(dt));
        pos = 0;
        for i=1:length(dt)
            if dt(i)>=delay
                size_steps=floor(tmax/dt(i)+offset);
            else
                size_steps=floor(tmax/delay+1+offset);
            end
            pos = pos + size_steps + 1;
            data(1:size_steps,1:3,i) = rawData(pos-(size_steps+1)+2:pos, 1:3);
        end

        step_sorted_data = [dt', reshape(data(end,2:3,:),[2,length(dt)])'];
        [Q(j), dQ(j)] = textra(step_sorted_data, 0, 0.1);

    end

    if(unitsFlag==1)
%         data(:,1,:) = data(:,1,:)/(4*sigma^2);
        %for Viscosity data
%         data(:,2:3,:) = data(:,2:3,:)/(4*sigma^2);
        Q = Q/(4*sigma^2);
        dQ = dQ/(4*sigma^2);
    elseif (unitsFlag==2)
%         data(:,1,:) = data(:,1,:)/(4*sigma^2);
        %for Psi1, Psi2 data
%         data(:,2:3,:) = data(:,2:3,:)/(16*sigma^4);
        Q = Q/(16*sigma^4);
        dQ = dQ/(16*sigma^4);
    elseif (unitsFlag==3)
%         data(:,1,:) = data(:,1,:)/(4*sigma^2);
    elseif (unitsFlag==4)
%         data(:,1,:) = data(:,1,:)/(4*sigma^2);
%         data(:,2:3,:) = abs(data(:,2:3,:));
        Q = abs(Q);
        dQ = abs(dQ);
    end
    
    if (plotFlag==1)
        figure();
        hold on
        axes1 = gca;
        axes1.XScale='log';
        axes1.YScale='log';
        errorbar(xvals, Q, dQ,...
                'DisplayName','TEXTRA','LineWidth',2);
        [~,~,~,~]=legend({},'Location','northwest',...
            'FontSize',16,'Interpreter','latex','Box','off');
    %     str = ['H = ', num2str(sigma^2); 'sr = ', num2str(];
    %     dim = [0.18 0.18 0.7 0.7];
    %     annotation('textbox',dim,'String',str,'FitBoxToText',...
    %         'on','Interpreter','latex','FontSize',14,'EdgeColor','None');
        title(regexprep(outputFileName,'\/.+\/',''))
        hold off
    end

    % If you want to save collected data to file
%     outputFileName;
    saveFile = fopen(outputFileName, 'w');
    fprintf(saveFile, '%g %g %g \n', [xvals; Q; dQ]);
    fclose(saveFile);
    
    f = [xvals; Q; dQ];
    
end
