function f = read_test_no_plot(inputFileName, outputFileName)

    dtData = dlmread('timestepdata.inp', '');
    dt = dtData(2:end, 1)';
    Nrelax_times = dtData(1, 3);
    delay = dtData(1,4);
    size_steps=floor(Nrelax_times/delay+1);
    inpdata = dlmread('inputparameters.inp', '');
    sigma = inpdata(4);
    alpha = inpdata(2);
    %set Rodlike_opt to true for Rodlike units, false for Hookean
    Rodlike_opt = true;

    rawData = dlmread(inputFileName, '');

    data = nan(size_steps,3,length(dt));
    pos = 0;
    for i=1:length(dt)
        if dt(i)>=delay
            size_steps=floor(Nrelax_times/dt(i));
        else
            size_steps=floor(Nrelax_times/delay+1);
        end
        pos = pos + size_steps + 1;
        data(1:size_steps,1:3,i) = rawData(pos-(size_steps+1)+2:pos, 1:3);
    end

    for i=1:size_steps
        step_sorted_data = [dt', reshape(data(i,2:3,:),[2,length(dt)])'];
        [Q(i), dQ(i)] = textra(step_sorted_data, 0, 0.1);
    end

    if(Rodlike_opt)
        data(:,1,:) = data(:,1,:)/(4*sigma^2);
        %for Viscosity data
        data(:,2:3,:) = data(:,2:3,:)/(4*sigma^2);
        Q = Q/(4*sigma^2);
        dQ = dQ/(4*sigma^2);
        %for Psi1, Psi2 data
    %     data(:,2:3,:) = data(:,2:3,:)/(16*sigma^4);
    %     Q = Q/(16*sigma^4);
    %     dQ = dQ/(16*sigma^4);

    end

    % If you want to save collected data to file
    outputFileName
    saveFile = fopen(outputFileName, 'w')
    fprintf(saveFile, '%g %g %g \n', [data(:,1,end)'; Q; dQ]);
end
