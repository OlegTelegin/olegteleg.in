% 1) This is the code to generate observations from https://olegtelegin.github.io/files/OlegTeleginQuietPeriod.pdf
% 2) Code is more straightforward to comprehend with the 
% code flowcharts: https://olegtelegin.github.io/files/OlegTeleginCodeFlowchartQuiet.pdf
% 3) It consists of 2 ingredients: ML-like algorithm to find a solution for
% one set of parameter values (1 observation), and Monte Carlo to accumulate
% observations. So, the code only accumulates observations without any statistical
% analysis due to the computational complexity of accumulating enough observations in 1 run. 
% 4) You can stop the code before n_outerrun (number of observations) is reached, generating fewer observations;
% then run the last line alone to write the results into .xlsx file (unfinished
% observation won't be written). Alternatively, you can set n_outerrun to a smaller number.
% 5) Displayed is the observation number in making (current n_outerrun).
% 6) All the comments, labeling the variables, use LaTeX mathematical
% symbols (like \beta, \sigma, ...), as in the paper.
% 7) The code is not optimized in the sense that it stores a lot of
% info in memory, not used in the final version but used to write/check
% the working code. However, this doesn't affect the speed by any significant amount.
% 8) Time spent to calculate 1 observation differs by a lot, but on average
% a mediocre laptop with i7proc/16GB RAM can calculate around  ??? observations per hour.

% number of draws (observations) 
n_outerrun = 1000;
% number of pairs i/d (\varepsilon_2 and \beta) for 1 draw to obtain
% characteristics of the risk premium (large model)
n_runs=10000;
% create empty table and set the size for the variables
table4 = table();
count_allruns_store = zeros(n_outerrun, 1);

% outer loop - for 1 draw
for outerrun = 1:n_outerrun
    % shuffle random number generator
    rng('shuffle');
    % set the size for the variables
    i = zeros(n_runs,1);
    d = zeros(n_runs,1);
    classres = zeros(n_runs,1);
    count1 = zeros(n_runs,1);
    count2 = zeros(n_runs,1);
    sum1ch = zeros(n_runs, 1);
    sum2ch = zeros(n_runs, 1);
    class_dop = ones(n_runs, 1);
    sum1_old2 = 0;
    % renew the table used within the loop
    second_table = table();
    % draw the parameters of the model
    a = 0.005; % CARA parameter
    b = rand(); % auxiliary variable for \delta_1, \delta_2, \delta_3
    c = rand(); % auxiliary variable for \delta_1, \delta_2, \delta_3
    e = rand(); % auxiliary variable for \sigma
    myArray = [0.01, 0.1, 1, 10, 100]; % o_1 possible values
    myArray1 = [0, 0.01, 0.1, 0.5, 1, 2, 10, 100]; % o_2 possible values
    randomIndex = randi(numel(myArray)); % pick the number from myArray
    randomIndex1 = randi(numel(myArray1)); % pick the number from my Array1
    o = myArray(randomIndex); % pick o_1
    odop = myArray1(randomIndex1); % pick o_2
    f = -1 + rand(); % \rho
    g = 10*rand(); % auxiliary variable for \sigma
    k = min(b, c); % \delta_1
    l = max(b, c)-min(b, c); % \delta_2
    m = 1-max(b, c); % \delta_3
    n = sqrt(e - g*log(rand())); % \sigma
    h = randn()*sqrt(min(b, c)); % \varepsilon_1
    mui = f*h*sqrt(l/k); % mean of \varepsilon_2
    sigmai = sqrt((1-f^2)*l); % standard deviation of \varepsilon_2
    p_i = makedist('Normal', mui, sigmai); % normal distribution of \varepsilon_2 to create percentiles
    j = randn()*sqrt(m); % \varepsilon_3
    % calculate percentiles for \varepsilon_2
    igrid3 = icdf(p_i, linspace(0,1,10001));
    xedges3 = igrid3;
    xedges3(1)=2*xedges3(2)-xedges3(3);
    xedges3(10001)=2*xedges3(10000)-xedges3(9999);
    % to keep track of the number of the draws (one less completed)
    disp(outerrun);
    
    % to obtain the initial guess for the classes closer to the equilibrium we start with
    % the small model
    % number of runs for the small model (to obtain stability)
    n_outerruninit = 50;
    % number of pairs i/d (\varepsilon_2 and \beta) for 1 draw (usually 400 for the small model)
    % to obtain characteristics of the risk premium
    n_runsinit=400;
    % create tables and variables
    count_allruns_storeinit = zeros(n_outerruninit, 1);
    second_tableinit = table();
    sumudiffinit_values = zeros(n_outerruninit, 1);
    
    % outer loop for the small model 
    for outerruninit = 1:n_outerruninit
        % set the size for the variables
        rng('shuffle');
        iinit = zeros(n_runsinit,1);
        dinit = zeros(n_runsinit,1);
        classresinit = zeros(n_runsinit,1);
        count1init = zeros(n_runsinit,1);
        count2init = zeros(n_runsinit,1);
        sum1chinit = zeros(n_runsinit, 1);
        sum2chinit = zeros(n_runsinit, 1);
        class_dopinit = ones(n_runsinit, 1);
        sum1_old2init = 0;
        
        % generate initial classes 
        for runinit = 1:n_runsinit
            % in the small model we start from the random initial guess,
            % assigning classes randomly for each (\varepsilon_2;\beta) dot
            classresinit(runinit) = randi([1,2]);
        end
        
        % uniformly put sqrt(n_runsinit)*sqrt(n_runsinit) grid of 
        % (\varepsilon_2;\beta) observations on this grid
        igrid2init = icdf(p_i, linspace(0,1,sqrt(n_runsinit)*50));
        xedges2init = igrid2init;
        for rungridinit1 = 1:sqrt(n_runsinit)
            for rungridinit2 = 1:sqrt(n_runsinit)
                mid_iinit = xedges2init(rungridinit1*50-25);
                mid_dinit = (2/sqrt(n_runsinit))*rungridinit2-(1/sqrt(n_runsinit));
                igrid2initdop = rungridinit2*sqrt(n_runsinit) - sqrt(n_runsinit) + rungridinit1;
                iinit(igrid2initdop)=mid_iinit;
                dinit(igrid2initdop)=mid_dinit;
            end
        end
        
        % calcucalte initial difference between U(with) and U(without) for the random classresinit
        count_onesinit = 0;
        count_allrunsinit = 0;
        meplus = mean(iinit(classresinit == 2));
        meminus = mean(iinit(classresinit == 1));
        mbplus = mean(dinit(classresinit == 2));
        mbminus = mean(dinit(classresinit == 1));
        veminus = std(iinit(classresinit == 1))*std(iinit(classresinit == 1));
        pr = mean(classresinit)-1;
        % calculate mbplustwo
        mbplustwo = zeros(n_runsinit, 1);
        selectedIndices = cell(n_runsinit, 1);
        for idop = 1:n_runsinit
            valuei = iinit(idop);
            [~, indices] = mink(abs(iinit - valuei), 3*sqrt(n_runsinit));
            selectedIndices{idop} = indices;
            selectedclass = classresinit(indices);
            selectedd = dinit(indices);
            mbplustwodop = mean(selectedd(selectedclass == 2));
            if isnan(mbplustwodop)
                mbplustwodop = 0;
            end
            mbplustwo(idop) = mbplustwodop;
        end
        % risk premium
        rp1 = a*(n^2)*(pr*mbplus*m+(1-pr)*(veminus+mbminus*m))+a*(n^2)*pr*(1-pr)*((meplus-meminus).^2)...
            +0.5*(a^3)*(n^4)*pr*(1-pr)*((mbplus*m-veminus-mbminus*m).^2)...
            -1.5*(a^2)*(n^3)*pr*m*mean(dinit(classresinit == 2))*mean(iinit(classresinit == 2))...
            -1.5*(a^2)*(n^3)*(1-pr)*(meminus*veminus+m*mean(dinit(classresinit == 1))*mean(iinit(classresinit == 1)))...
            +1.5*(a^2)*(n^3)*(pr*meplus+(1-pr)*meminus)*(pr*mbplus*m+(1-pr)*(veminus+mbminus*m));
        % U(with) - U(without)
        udiffinit = o*((dinit.*(n^2)*m-(n^2)*mbplustwo.*m).^2)...
            +odop*((iinit.*n-a*(n^2)*mbplustwo.*m-n*pr*meplus-n*(1-pr)*meminus+rp1).^2)...
            +((a*(n^2)*mbplustwo.*m).^2)-o*((dinit.*(n^2)*m-(n^2)*veminus-(n^2)*mbminus*m).^2)...
            -odop*((a*(n^2)*(veminus+mbminus*m)+n*pr*meplus-n*pr*meminus-rp1).^2)...
            -((iinit.*n-n*meminus+a*(n^2)*(veminus+mbminus*m)).^2);
        
        % calculate mistakes (count1init) and count mistakes (count2init)
        classcorrectinit = zeros(n_runsinit,1);
        for runinit = 1:n_runsinit
            classcorrectinit(runinit) = (1.5-classresinit(runinit))*udiffinit(runinit);
            if classcorrectinit(runinit) >=0
                count1init(runinit) = 0;
                count2init(runinit) = 0;
            else
                count1init(runinit) = abs(udiffinit(runinit));
                count2init(runinit) = 1;
            end
        end
        
        % store classresinit/count1/count2 (further will be used as the values on
        % the previous run)
        class_startinit = classresinit;
        count1_startinit = count1init;
        count2_startinit = count2init;
        
        % calculate sum/number of mistakes
        sum1_oldinit = sum(count1init);
        sum2_oldinit = sum(count2init);
        
        % loop to reassign class to get rid of mistakingly assigned classes
        while true
            % find the index of observation with maximal mistake and change this index
            max_count1init = max(count1init);
            idx3init = find(count1init == max_count1init);
            classresinit(idx3init) = 3 - classresinit(idx3init);
            
            % calcucalte the difference between U(with) and U(without) again
            meplus = mean(iinit(classresinit == 2));
            meminus = mean(iinit(classresinit == 1));
            mbplus = mean(dinit(classresinit == 2));
            mbminus = mean(dinit(classresinit == 1));
            veminus = std(iinit(classresinit == 1))*std(iinit(classresinit == 1));
            if isnan(veminus)
                veminus = 0;
            end
            pr = mean(classresinit)-1;
            proizvplus = mean(dinit(classresinit == 2))*mean(iinit(classresinit == 2));
            proizvminus = mean(dinit(classresinit == 1))*mean(iinit(classresinit == 1));
            if isnan(meplus)
                meplus = 0;
                mbplus = 0;
                proizvplus = 0;
            end
            if isnan(meminus)
                meminus = 0;
                mbminus = 0;
                proizvminus = 0;
            end
            % calculate mbplustwo
            mbplustwo = zeros(n_runsinit, 1);
            for idop = 1:n_runsinit
                valuei = iinit(idop);
                indices = selectedIndices{idop};
                selectedclass = classresinit(indices);
                selectedd = dinit(indices);
                mbplustwodop = mean(selectedd(selectedclass == 2));
                if isnan(mbplustwodop)
                    mbplustwodop = 0;
                end
                mbplustwo(idop) = mbplustwodop;
            end
            % risk premium
            rp1 = a*(n^2)*(pr*mbplus*m+(1-pr)*(veminus+mbminus*m))+a*(n^2)*pr*(1-pr)*((meplus-meminus).^2)...
                +0.5*(a^3)*(n^4)*pr*(1-pr)*((mbplus*m-veminus-mbminus*m).^2)...
                -1.5*(a^2)*(n^3)*pr*m*proizvplus-1.5*(a^2)*(n^3)*(1-pr)*(meminus*veminus+m*proizvminus)...
                +1.5*(a^2)*(n^3)*(pr*meplus+(1-pr)*meminus)*(pr*mbplus*m+(1-pr)*(veminus+mbminus*m));
            % U(with) - U(without)
            udiffinit = o*((dinit.*(n^2)*m-(n^2)*mbplustwo.*m).^2)...
                +odop*((iinit.*n-a*(n^2)*mbplustwo.*m-n*pr*meplus-n*(1-pr)*meminus+rp1).^2)...
                +((a*(n^2)*mbplustwo.*m).^2)-o*((dinit.*(n^2)*m-(n^2)*veminus-(n^2)*mbminus*m).^2)...
                -odop*((a*(n^2)*(veminus+mbminus*m)+n*pr*meplus-n*pr*meminus-rp1).^2)...
                -((iinit.*n-n*meminus+a*(n^2)*(veminus+mbminus*m)).^2);
            
            % calculate mistakes (count1init) and count mistakes (count2init)
            classcorrectinit = zeros(n_runsinit,1);
            for runinit = 1:n_runsinit
                classcorrectinit(runinit) = (1.5-classresinit(runinit))*udiffinit(runinit);
                if classcorrectinit(runinit) >=0
                    count1init(runinit) = 0;
                    count2init(runinit) = 0;
                else
                    count1init(runinit) = abs(udiffinit(runinit));
                    count2init(runinit) = 1;
                end
            end
            
            % counter for count2init==1 - we will further escape the loop if we can't
            % progress to 0; and counter for the number of runs
            if (sum(count2init)==1)
                count_onesinit = count_onesinit+1;
            end
            count_allrunsinit = count_allrunsinit+1;
            count_allruns_storeinit(outerruninit) = count_allrunsinit;
            
            min_mistakeinit = min(sum(count2init), sum2_oldinit);
            
            % stop if no mistakes, or we loop on 1 mistake, or we hit the limit of 10000 runs
            % (in this case we are within the loop of the length n, say 17 12 5 3 9 17 ...
            % mistakes, so we stop on 3, also checking that this number of mistakes doesn't exceed
            % sqrt(n_runsinit)), or we hit the limit of 11000 runs (happens when
            % in a previous case we have a loop of say 4 4 4 4 4... mistakes)
            % practically for the smaller model count2init=0 all the time,
            % this block of conditions is needed only for the large model (n_runs=10000)
            if (sum(count2init) < 1) || ...
                    (count_onesinit == 10) || ...
                    (count_allrunsinit > 2999 && sum(count2init) > sum2_oldinit && ...
                    min_mistakeinit < sqrt(n_runsinit)) || ...
                    (count_allrunsinit == 4000)
                break;
            end
            
            % update sum1 and sum2 in the original table
            sum1_old2init = sum1_oldinit;
            sum1_oldinit = sum(count1init);
            sum2_oldinit = sum(count2init);
            class_startinit = classresinit;
            count1_startinit = count1init;
            count2_startinit = count2init;
            
        end
        
        % calculate the sum of udiff for all the observations
        absudiffinit = abs(udiffinit);
        sumudiffinit = sum(absudiffinit);
        sumudiffinit_values(outerruninit) = sumudiffinit;
        
        % store draws in a table    
        second_tableinit = [second_tableinit; table(iinit,dinit,classresinit, udiffinit, count1init, ...
            count2init, sum1chinit, sum2chinit, class_startinit, count1_startinit, count2_startinit)];
        
    end
    
    % calculate average classes from the small model
    reshapedMatrix = reshape(second_tableinit.classresinit, n_runsinit, n_outerruninit);
    numSmallest = 0.5*n_outerruninit;
    [sortedValues, sortedIndices] = sort(sumudiffinit_values, 'ascend');
    indicesOfSmallest = sortedIndices(1:numSmallest);
    reshapedMatrix1 = reshapedMatrix;
    reshapedMatrix1(:, indicesOfSmallest) = [];
    summedclassresinit = sum(reshapedMatrix1, 2);
    threshold = 0.75*n_outerruninit;
    avclassresinit = ones(size(summedclassresinit));
    avclassresinit(summedclassresinit >= threshold) = 2;
    
    % the main loop is back - we draw (i;d) for the large model, now randomly, not on an even grid
    for run = 1:n_runs
        i(run) = random(p_i);
        d(run) = 2*rand();
    end
    
    % now we draw the initial guess for classres according to the small
    % model solution - divide (\varepsilon_2;\beta) space into an equal grid, where the
    % observations from the small model are in the center of the grid cells;
    % then assign the class from the small model to all the dots in the grid cell
    binEdges = icdf(p_i, linspace(0,1,sqrt(n_runsinit)+1));
    binEdges(1) = min(i)-10;
    binEdges(sqrt(n_runsinit)+1) = max(i)+10;
    binEdgesvert = linspace(0,2,sqrt(n_runsinit)+1);
    horindices = zeros(n_runs,1);
    vertindices = zeros(n_runs,1);
    
    % assign the value of the classresinit to the whole grid cell
    for run = 1:n_runs
        horindices(run) = find(i(run) <= binEdges, 1)-1;
        vertindices(run) = find(d(run) <= binEdgesvert, 1)-1;
        final_index = vertindices(run)*sqrt(n_runsinit)-sqrt(n_runsinit)+horindices(run);
        classres(run) = avclassresinit(final_index);
    end
    
    % calcucalte the difference between U(with) and U(without)
    count_ones = 0;
    count_allruns = 0;
    meplus = mean(i(classres == 2));
    meminus = mean(i(classres == 1));
    mbplus = mean(d(classres == 2));
    mbminus = mean(d(classres == 1));
    veminus = std(i(classres == 1))*std(i(classres == 1));
    pr = mean(classres)-1;
    % calculate mbplustwo
    mbplustwo = zeros(n_runs, 1);
    selectedIndices = cell(n_runs, 1);
    for idop = 1:n_runs
        valuei = i(idop);
        [~, indices] = mink(abs(i - valuei), 200);
        selectedIndices{idop} = indices;
        selectedclass = classres(indices);
        selectedd = d(indices);
        mbplustwodop = mean(selectedd(selectedclass == 2));
        if isnan(mbplustwodop)
            mbplustwodop = 0;
        end
        mbplustwo(idop) = mbplustwodop;
    end
    % risk premium
    rp1 = a*(n^2)*(pr*mbplus*m+(1-pr)*(veminus+mbminus*m))+a*(n^2)*pr*(1-pr)*((meplus-meminus).^2)...
        +0.5*(a^3)*(n^4)*pr*(1-pr)*((mbplus*m-veminus-mbminus*m).^2)...
        -1.5*(a^2)*(n^3)*pr*m*mean(d(classres == 2))*mean(i(classres == 2))...
        -1.5*(a^2)*(n^3)*(1-pr)*(meminus*veminus+m*mean(d(classres == 1))*mean(i(classres == 1)))...
        +1.5*(a^2)*(n^3)*(pr*meplus+(1-pr)*meminus)*(pr*mbplus*m+(1-pr)*(veminus+mbminus*m));
    % U(with) - U(without)
    udiff = o*((d.*(n^2)*m-(n^2)*mbplustwo.*m).^2)...
        +odop*((i.*n-a*(n^2)*mbplustwo.*m-n*pr*meplus-n*(1-pr)*meminus+rp1).^2)...
        +((a*(n^2)*mbplustwo.*m).^2)-o*((d.*(n^2)*m-(n^2)*veminus-(n^2)*mbminus*m).^2)...
        -odop*((a*(n^2)*(veminus+mbminus*m)+n*pr*meplus-n*pr*meminus-rp1).^2)...
        -((i.*n-n*meminus+a*(n^2)*(veminus+mbminus*m)).^2);
    
    % calculate mistakes (count1) and count mistakes (count2)
    classcorrect = zeros(n_runs,1);
    for run = 1:n_runs
        classcorrect(run) = (1.5-classres(run))*udiff(run);
        if classcorrect(run) >=0
            count1(run) = 0;
            count2(run) = 0;
        else
            count1(run) = abs(udiff(run));
            count2(run) = 1;
        end
    end
    
    % store classres/count1/count2 (further will be used as the values on the previous run)
    class_start = classres;
    count1_start = count1;
    count2_start = count2;
    sum1_old = sum(count1);
    sum2_old = sum(count2);
    
    % loop to reassign class to get rid of mistakingly assigned classes
    while true
        % find the index of observation with maximal mistake and change this index
        max_count1 = max(count1);
        idx3 = find(count1 == max_count1);
        classres(idx3) = 3 - classres(idx3);
        
        % calcucalte U(with) - U(without) with new classes
        meplus = mean(i(classres == 2));
        meminus = mean(i(classres == 1));
        mbplus = mean(d(classres == 2));
        mbminus = mean(d(classres == 1));
        veminus = std(i(classres == 1))*std(i(classres == 1));
        if isnan(veminus)
            veminus = 0;
        end
        pr = mean(classres)-1;
        proizvplus = mean(d(classres == 2))*mean(i(classres == 2));
        if isnan(meplus)
            meplus = 0;
            mbplus = 0;
            proizvplus = 0;
        end
        proizvminus = mean(d(classres == 1))*mean(i(classres == 1));
        if isnan(meminus)
            meminus = 0;
            mbminus = 0;
            proizvminus = 0;
        end
        % calculate mbplustwo
        mbplustwo = zeros(n_runs, 1);
        for idop = 1:n_runs
            valuei = i(idop);
            indices = selectedIndices{idop};
            selectedclass = classres(indices);
            selectedd = d(indices);
            mbplustwodop = mean(selectedd(selectedclass == 2));
            if isnan(mbplustwodop)
                mbplustwodop = 0;
            end
            mbplustwo(idop) = mbplustwodop;
        end
        % risk premium
        rp1 = a*(n^2)*(pr*mbplus*m+(1-pr)*(veminus+mbminus*m))+a*(n^2)*pr*(1-pr)*((meplus-meminus).^2)...
            +0.5*(a^3)*(n^4)*pr*(1-pr)*((mbplus*m-veminus-mbminus*m).^2)...
            -1.5*(a^2)*(n^3)*pr*m*proizvplus-1.5*(a^2)*(n^3)*(1-pr)*(meminus*veminus+m*proizvminus)...
            +1.5*(a^2)*(n^3)*(pr*meplus+(1-pr)*meminus)*(pr*mbplus*m+(1-pr)*(veminus+mbminus*m));
        % U(with) - U(without)
        udiff = o*((d.*(n^2)*m-(n^2)*mbplustwo.*m).^2)...
            +odop*((i.*n-a*(n^2)*mbplustwo.*m-n*pr*meplus-n*(1-pr)*meminus+rp1).^2)...
            +((a*(n^2)*mbplustwo.*m).^2)-o*((d.*(n^2)*m-(n^2)*veminus-(n^2)*mbminus*m).^2)...
            -odop*((a*(n^2)*(veminus+mbminus*m)+n*pr*meplus-n*pr*meminus-rp1).^2)...
            -((i.*n-n*meminus+a*(n^2)*(veminus+mbminus*m)).^2);
        
        % calculate mistakes (count1) and count mistakes (count2)
        classcorrect = zeros(n_runs,1);
        for run = 1:n_runs
            classcorrect(run) = (1.5-classres(run))*udiff(run);
            if classcorrect(run) >=0
                count1(run) = 0;
                count2(run) = 0;
            else
                count1(run) = abs(udiff(run));
                count2(run) = 1;
            end
        end
        
        % counter for count2init==1 - we will further escape the loop if we can't
        % progress to 0; and counter for the number of runs
        if (sum(count2)==1)
            count_ones = count_ones+1;
        end
        count_allruns = count_allruns+1;
        count_allruns_store(outerrun) = count_allruns;
        
        min_mistake = min(sum(count2), sum2_old);
        
        % stop if no mistakes, or we looped on 1 mistake, or we hit the limit of 10000 runs
        % (in this case we are within the loop, say 17 12 5 3 9 17 ...
        % mistakes, so we stop on 3, also checking that this number of mistakes doesn't exceed
        % 30), or we hit the limit of 11000 runs (happens when
        % in a previous case we have a loop of say 4 4 4 4 4... mistakes)
        if (sum(count2)<1) || (count_ones == 10) || ...
                (count_allruns > 9999 && sum(count2)>sum2_old && min_mistake<30) || ...
                (count_allruns == 11000)
            break;
        end
        
        % update sum1 and sum2 in the original table
        sum1_old2 = sum1_old;
        sum1_old = sum(count1);
        sum2_old = sum(count2);
        class_start = classres;
        count1_start = count1;
        count2_start = count2;
        
    end
    
    % we decide which classes to take - from the last run or from the second
    % to last, and then we store the results
    if (sum(count2)>sum2_old)
        classtotake = class_start;
    else
        classtotake = classres;
    end

    % store draws in a table    
    second_table = [second_table; table(i,d,classres, udiff, count1, count2, sum1ch, sum2ch, ...
        class_start, count1_start, count2_start, mbplustwo, classtotake)];
    
    store_eminus = mean(second_table.i(second_table.classtotake==1));
    store_eplus = mean(second_table.i(second_table.classtotake==2));
    store_bminus = mean(second_table.d(second_table.classtotake==1));
    store_bplus = mean(second_table.d(second_table.classtotake==2));
    store_pr = mean(second_table.classtotake)-1;
    store_veminus = std(second_table.i(second_table.classtotake==1))*std(second_table.i(second_table.classtotake==1));
    store_vbminus = std(second_table.d(second_table.classtotake==1))*std(second_table.d(second_table.classtotake==1));
    store_vbplus = std(second_table.d(second_table.classtotake==2))*std(second_table.d(second_table.classtotake==2));
    store_proizvminus = mean(second_table.d(second_table.classtotake==1))*mean(second_table.i(second_table.classtotake==1));
    store_proizvplus = mean(second_table.d(second_table.classtotake==2))*mean(second_table.i(second_table.classtotake==2)); 
    store_bplustwo_dop2 = zeros(10000, 1);
    for rungraph7 = 1:10000
        if (rungraph7>20)
            min_bound = rungraph7-20;
        else
            min_bound = 1;
        end
        if (rungraph7<9980)
            max_bound = rungraph7+21;
        else
            max_bound = 10000;
        end
        store_bplustwo_dop = mean(second_table.d(second_table.classtotake==2 & ...
            second_table.i>xedges3(min_bound) & second_table.i<xedges3(max_bound)));
        if isnan(store_bplustwo_dop)
            store_bplustwo_dop = 0;
        end
        store_bplustwo_dop2(rungraph7) = store_bplustwo_dop;
    end
    store_bplustwo = smoothdata(store_bplustwo_dop2, 'gaussian', 1000);
    store_rp1 = a*(n^2)*(store_pr*store_bplus*m...
        +(1-store_pr)*(store_veminus+m*store_bminus))...
        +a*(n^2)*store_pr*(1-store_pr)*((store_eplus-store_eminus)^2)...
        +0.5*(a^3)*(n^4)*store_pr*(1-store_pr)*((store_bplus*m-store_veminus-store_bminus*m)^2)...
        -1.5*(a^2)*(n^3)*store_pr*m*store_proizvplus...
        -1.5*(a^2)*(n^3)*(1-store_pr)*(store_eminus*store_veminus+m*store_proizvminus)...
        +1.5*(a^2)*(n^3)*(store_pr*store_eplus...
        +(1-store_pr)*store_eminus)*(store_pr*store_bplus*m...
        +(1-store_pr)*(store_veminus+store_bminus*m));
    
    % set the 10000*10000 grid for the second step of the algorithm
    dgrid2 = linspace(0,2,10001);
    igrid2 = icdf(p_i, linspace(0,1,10001));
    [X2,Y2] = meshgrid(igrid2,dgrid2);
    xedges2 = igrid2;
    xedges2(1)=2*xedges2(2)-xedges2(3);
    xedges2(10001)=2*xedges2(10000)-xedges2(9999);
    yedges2 = dgrid2;  
    
    % replace omitted values for the case store_pr == 1
    if (store_pr == 1)
        store_eminus = 0;
        store_bminus = 0;
        store_veminus = 0;
        store_vbminus = 0;
        store_proizvminus = 0;
        store_rp1 = a*(n^2)*(store_pr*store_bplus*m+(1-store_pr)*(store_veminus+store_bminus*m))...
            +a*(n^2)*store_pr*(1-store_pr)*((store_eplus-store_eminus)^2)...
            +0.5*(a^3)*(n^4)*store_pr*(1-store_pr)*((store_bplus*m-store_veminus-store_bminus*m)^2)...
            -1.5*(a^2)*(n^3)*store_pr*m*store_proizvplus...
            -1.5*(a^2)*(n^3)*(1-store_pr)*(store_eminus*store_veminus+m*store_proizvminus)...
            +1.5*(a^2)*(n^3)*(store_pr*store_eplus+(1-store_pr)*store_eminus)*(store_pr*store_bplus*m...
            +(1-store_pr)*(store_veminus+store_bminus*m));
    end
    
    % replace omitted values for the case store_pr == 0
    if (store_pr == 0)
        store_eplus = 0;
        store_bplus = 0;
        store_vbplus = 0;
        store_proizvplus = 0;
        store_rp1 = a*(n^2)*(store_pr*store_bplus*m+(1-store_pr)*(store_veminus+store_bminus*m))...
            +a*(n^2)*store_pr*(1-store_pr)*((store_eplus-store_eminus)^2)...
            +0.5*(a^3)*(n^4)*store_pr*(1-store_pr)*((store_bplus*m-store_veminus-store_bminus*m)^2)...
            -1.5*(a^2)*(n^3)*store_pr*m*store_proizvplus...
            -1.5*(a^2)*(n^3)*(1-store_pr)*(store_eminus*store_veminus+m*store_proizvminus)...
            +1.5*(a^2)*(n^3)*(store_pr*store_eplus...
            +(1-store_pr)*store_eminus)*(store_pr*store_bplus*m+(1-store_pr)*(store_veminus+store_bminus*m));
    end
    
    % calcucalte U(with) - U(without) for each dot and store 0 or 1 for the class and Udiff itself
    final_grid = zeros(length(igrid2)-1,length(dgrid2)-1);
    final_grid2 = zeros(length(igrid2)-1,length(dgrid2)-1);
    for rungraph3 = 1:size(X2,1)-1
        for rungraph4 = 1:size(X2,2)-1
            mid_i = (xedges2(rungraph3+1)+xedges2(rungraph3))/2;
            mid_d = (yedges2(rungraph4+1)+yedges2(rungraph4))/2;
            store_bplustwo_formula = store_bplustwo(rungraph3);
            udiff_final = o*((mid_d*(n^2)*m-(n^2)*store_bplustwo_formula*m)^2)...
                +odop*((mid_i*n-a*(n^2)*store_bplustwo_formula*m-n*store_pr*store_eplus-n*(1-store_pr)*store_eminus+store_rp1)^2)...
                +((a*(n^2)*store_bplustwo_formula*m)^2)-o*((mid_d*(n^2)*m-(n^2)*store_veminus-(n^2)*store_bminus*m)^2)...
                -odop*((a*(n^2)*(store_veminus+store_bminus*m)+n*store_pr*store_eplus-n*store_pr*store_eminus-store_rp1)^2)...
                -((mid_i*n-n*store_eminus+a*(n^2)*(store_veminus+store_bminus*m))^2);
            if udiff_final>0
                udiff_finsign=0;
            else
                udiff_finsign=1;
            end
            final_grid(rungraph3,rungraph4) = udiff_finsign;
            final_grid2(rungraph3,rungraph4) = udiff_final;
        end
    end
    
    % now we store all the equilibrium characteristics of the observation
    mean_final_grid2 = mean(final_grid2(:));
    var_final_grid2 = var(final_grid2(:));
    num_ones = sum(final_grid(:) == 1);
    if (count_allruns_store(outerrun) > 9999)
        infin_store=1;
    else
        infin_store=0;
    end
    proportion_ones = num_ones / numel(final_grid);
    
    for i = 1:10
        column_indices = i * 1000;
        eval(['num_ones' num2str(i) ' = sum(final_grid(:, column_indices) == 1);']);
    end
    num_ones0 = sum(final_grid(:,1) == 1);
    
    propminustwosigma = sum(final_grid(228,:) == 1) / size(final_grid,2);
    propminussigma = sum(final_grid(1587,:) == 1) / size(final_grid,2);
    propminushalfsigma = sum(final_grid(3086,:) == 1) / size(final_grid,2);
    propminusquartersigma = sum(final_grid(4013,:) == 1) / size(final_grid,2);
    propminuszerosigma = (sum(final_grid(5000,:) == 1)+sum(final_grid(5001,:) == 1))/(2*size(final_grid,2));
    propplusquartersigma = sum(final_grid(5988,:) == 1) / size(final_grid,2);
    propplushalfsigma = sum(final_grid(6915,:) == 1) / size(final_grid,2);
    propplussigma = sum(final_grid(8414,:) == 1) / size(final_grid,2);
    propplustwosigma = sum(final_grid(9773,:) == 1) / size(final_grid,2);
    percentiles_eminus = 10000*normcdf(store_eminus, mui, sigmai);
    percentiles_eplus = 10000*normcdf(store_eplus, mui, sigmai);
    
    if (num_ones == 100000000)
        % store all the characteristics for the case Prop^+==1
        zero_indices = NaN;
        one_indices = find(final_grid == 1);
        min_i = NaN;
        max_i = NaN;
        min_d = NaN;
        max_d = NaN;
        min_i1 = min(mod(one_indices - 1, size(final_grid, 1)) + 1);
        max_i1 = max(mod(one_indices - 1, size(final_grid, 1)) + 1);
        min_d1 = min(floor((one_indices - 1) / size(final_grid, 1)) + 1);
        max_d1 = max(floor((one_indices - 1) / size(final_grid, 1)) + 1);
        row_indices = NaN;
        avg_row = NaN;
        col_indices = NaN;
        avg_col = NaN;
        row_indices1 = floor((one_indices - 1) / size(final_grid, 1)) + 1;
        avg_row1 = mean(row_indices1);
        col_indices1 = mod(one_indices - 1, size(final_grid, 1)) + 1;
        avg_col1 = mean(col_indices1);
        variance_d_indices = NaN;
        variance_d_indices1 = var(row_indices1);
        variance_i_indices = NaN;
        variance_i_indices1 = var(col_indices1);
        % add another line to the table with all the draws
        table4 = [table4; table(proportion_ones, min_i, max_i, min_d, max_d, min_i1, max_i1, ...
            min_d1, max_d1, variance_d_indices, variance_d_indices1, variance_i_indices, ...
            variance_i_indices1, num_ones10, num_ones9, ...
            num_ones8, num_ones7, num_ones6, num_ones5, num_ones4, ...
            num_ones3, num_ones2, num_ones1, num_ones0, avg_row, avg_col, avg_row1, avg_col1, ...
            mean_final_grid2, var_final_grid2, propminustwosigma, propminussigma, ...
            propminushalfsigma, propminusquartersigma, propminuszerosigma, propplusquartersigma, ...
            propplushalfsigma, propplussigma, propplustwosigma, a, e, f, g, h, j, k, l, m, n, o, ...
            odop, mui, sigmai, store_eminus, store_eplus, store_bminus, store_bplus, store_pr, ...
            store_veminus, store_vbminus, store_vbplus, store_proizvminus, store_proizvplus, ...
            store_rp1, percentiles_eminus, percentiles_eplus, infin_store, min_mistake, count_allruns)];
    elseif (num_ones == 0)
        % store all the equilibrium characteristics for the case Prop^+==0
        zero_indices = find(final_grid == 0);
        one_indices = NaN;
        min_i = min(mod(zero_indices - 1, size(final_grid, 1)) + 1);
        max_i = max(mod(zero_indices - 1, size(final_grid, 1)) + 1);
        min_d = min(floor((zero_indices - 1) / size(final_grid, 1)) + 1);
        max_d = max(floor((zero_indices - 1) / size(final_grid, 1)) + 1);
        min_i1 = NaN;
        max_i1 = NaN;
        min_d1 = NaN;
        max_d1 = NaN;
        row_indices = floor((zero_indices - 1) / size(final_grid, 1)) + 1;
        avg_row = mean(row_indices);
        col_indices = mod(zero_indices - 1, size(final_grid, 1)) + 1;
        avg_col = mean(col_indices);
        row_indices1 = NaN;
        avg_row1 = NaN;
        col_indices1 = NaN;
        avg_col1 = NaN;
        variance_d_indices = var(row_indices);
        variance_d_indices1 = NaN;
        variance_i_indices = var(col_indices);
        variance_i_indices1 = NaN;
        % add another line to the table with all the draws
        table4 = [table4; table(proportion_ones, min_i, max_i, min_d, max_d, min_i1, max_i1, ...
            min_d1, max_d1, variance_d_indices, variance_d_indices1, variance_i_indices, ...
            variance_i_indices1, num_ones10, num_ones9, num_ones8, ...
            num_ones7, num_ones6, num_ones5, num_ones4, num_ones3, ...
            num_ones2, num_ones1, num_ones0, avg_row, avg_col, avg_row1, avg_col1, ...
            mean_final_grid2, var_final_grid2, propminustwosigma, propminussigma, ...
            propminushalfsigma, propminusquartersigma, propminuszerosigma, ...
            propplusquartersigma, propplushalfsigma, propplussigma, propplustwosigma, ...
            a, e, f, g, h, j, k, l, m, n, o, odop, mui, sigmai, store_eminus, ...
            store_eplus, store_bminus, store_bplus, store_pr, store_veminus, ...
            store_vbminus, store_vbplus, store_proizvminus, store_proizvplus, ...
            store_rp1, percentiles_eminus, percentiles_eplus, infin_store, min_mistake, count_allruns)];
    else
        % store all the equilibrium characteristics for Prop^+ != 0;1
        zero_indices = find(final_grid == 0);
        one_indices = find(final_grid == 1);
        min_i = min(mod(zero_indices - 1, size(final_grid, 1)) + 1);
        max_i = max(mod(zero_indices - 1, size(final_grid, 1)) + 1);
        min_d = min(floor((zero_indices - 1) / size(final_grid, 1)) + 1);
        max_d = max(floor((zero_indices - 1) / size(final_grid, 1)) + 1);
        min_i1 = min(mod(one_indices - 1, size(final_grid, 1)) + 1);
        max_i1 = max(mod(one_indices - 1, size(final_grid, 1)) + 1);
        min_d1 = min(floor((one_indices - 1) / size(final_grid, 1)) + 1);
        max_d1 = max(floor((one_indices - 1) / size(final_grid, 1)) + 1);
        row_indices = floor((zero_indices - 1) / size(final_grid, 1)) + 1;
        avg_row = mean(row_indices);
        col_indices = mod(zero_indices - 1, size(final_grid, 1)) + 1;
        avg_col = mean(col_indices);
        row_indices1 = floor((one_indices - 1) / size(final_grid, 1)) + 1;
        avg_row1 = mean(row_indices1);
        col_indices1 = mod(one_indices - 1, size(final_grid, 1)) + 1;
        avg_col1 = mean(col_indices1);
        variance_d_indices = var(row_indices);
        variance_d_indices1 = var(row_indices1);
        variance_i_indices = var(col_indices);
        variance_i_indices1 = var(col_indices1);
        % add another line to the table with all the draws
        table4 = [table4; table(proportion_ones, min_i, max_i, min_d, max_d, min_i1, ...
            max_i1, min_d1, max_d1, variance_d_indices, variance_d_indices1, ...
            variance_i_indices, variance_i_indices1, num_ones10, num_ones9, num_ones8, ...
            num_ones7, num_ones6, num_ones5, num_ones4, num_ones3, num_ones2, ...
            num_ones1, num_ones0, avg_row, avg_col, avg_row1, avg_col1, ...
            mean_final_grid2, var_final_grid2, propminustwosigma, propminussigma, ...
            propminushalfsigma, propminusquartersigma, propminuszerosigma, ...
            propplusquartersigma, propplushalfsigma, propplussigma, propplustwosigma, ...
            a, e, f, g, h, j, k, l, m, n, o, odop, mui, sigmai, store_eminus, ...
            store_eplus, store_bminus, store_bplus, store_pr, store_veminus, ...
            store_vbminus, store_vbplus, store_proizvminus, store_proizvplus, ...
            store_rp1, percentiles_eminus, percentiles_eplus, infin_store, min_mistake, count_allruns)];
    end
end

% write the table into .xlsx file
writetable(table4, 'simulation_results.xlsx');