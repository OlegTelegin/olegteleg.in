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
% a mediocre laptop with i7proc/16GB RAM can calculate around  15-20 observations per hour.

% number of draws (observations) 
n_outerrun = 1000;
% number of pairs i/d (\varepsilon_2 and \beta) for 1 draw to obtain
% characteristics of the risk premium (large model)
n_runs=10000;
% create empty table and set the size for the variables
table4 = table();
check_no = zeros(n_outerrun, 1);
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
    a_values = zeros(n_runs, 1);
    b_values = zeros(n_runs, 1);
    c_values = zeros(n_runs, 1);
    d_values = zeros(n_runs, 1);
    e_values = zeros(n_runs, 1);
    f_values = zeros(n_runs, 1);
    g_values = zeros(n_runs, 1);
    h_values = zeros(n_runs, 1);
    i_values = zeros(n_runs, 1);
    j_values = zeros(n_runs, 1);
    k_values = zeros(n_runs, 1);
    l_values = zeros(n_runs, 1);
    m_values = zeros(n_runs, 1);
    n_values = zeros(n_runs, 1);
    o_values = zeros(n_runs, 1);
    odop_values = zeros(n_runs, 1);
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
    check_noinit = zeros(n_outerruninit, 1);
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
        a_valuesinit = zeros(n_runsinit, 1);
        b_valuesinit = zeros(n_runsinit, 1);
        c_valuesinit = zeros(n_runsinit, 1);
        d_valuesinit = zeros(n_runsinit, 1);
        e_valuesinit = zeros(n_runsinit, 1);
        f_valuesinit = zeros(n_runsinit, 1);
        g_valuesinit = zeros(n_runsinit, 1);
        h_valuesinit = zeros(n_runsinit, 1);
        i_valuesinit = zeros(n_runsinit, 1);
        j_valuesinit = zeros(n_runsinit, 1);
        k_valuesinit = zeros(n_runsinit, 1);
        l_valuesinit = zeros(n_runsinit, 1);
        m_valuesinit = zeros(n_runsinit, 1);
        n_valuesinit = zeros(n_runsinit, 1);
        o_valuesinit = zeros(n_runsinit, 1);
        odop_valuesinit = zeros(n_runsinit, 1);
        sum1chinit = zeros(n_runsinit, 1);
        sum2chinit = zeros(n_runsinit, 1);
        class_dopinit = ones(n_runsinit, 1);
        sum1_old2init = 0;
        
        % generate initial classes 
        for runinit = 1:n_runsinit
            % in the small model we start from the random initial guess,
            % assigning classes randomly for each (\varepsilon_2;\beta) dot
            classresinit(runinit) = randi([1,2]);
            a_valuesinit(runinit) = a;
            e_valuesinit(runinit) = e;
            f_valuesinit(runinit) = f;
            g_valuesinit(runinit) = g;
            h_valuesinit(runinit) = h;
            j_valuesinit(runinit) = j;
            k_valuesinit(runinit) = k;
            l_valuesinit(runinit) = l;
            m_valuesinit(runinit) = m;
            n_valuesinit(runinit) = n;
            o_valuesinit(runinit) = o;
            odop_valuesinit(runinit) = odop;
        end
        
        % set the horizontal grid with equal percentiles
        igrid2init = icdf(p_i, linspace(0,1,sqrt(n_runsinit)*50));
        xedges2init = igrid2init;
        
        % uniformly put sqrt(n_runsinit)*sqrt(n_runsinit) grid of 
        % (\varepsilon_2;\beta) observations on this grid
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
        vbplus = std(dinit(classresinit == 2))*std(dinit(classresinit == 2));
        vbminus = std(dinit(classresinit == 1))*std(dinit(classresinit == 1));
        pr = mean(classresinit)-1;
        % calculate mbplustwo
        mbplustwo = zeros(n_runsinit, 1);
        selectedIndices = cell(n_runsinit, 1);
        % iterate through each row in the table
        for idop = 1:n_runsinit
            % get the value from column in the current row
            valuei = iinit(idop);
            % find the indices of the X closest values to the current value in column
            [~, indices] = mink(abs(iinit - valuei), 3*sqrt(n_runsinit));
            selectedIndices{idop} = indices;
            selectedclass = classresinit(indices);
            selectedd = dinit(indices);
            % count the number of 2s in column for the selected indices
            mbplustwodop = mean(selectedd(selectedclass == 2));
            if isnan(mbplustwodop)
                mbplustwodop = 0;
            end
            % assign the count to column in the current row
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
        for runinit = 1:n_runsinit
            if classresinit(runinit) == 1
                if udiffinit(runinit) >= 0
                    count1init(runinit) = 0;
                    count2init(runinit) = 0;
                else
                    count1init(runinit) = abs(udiffinit(runinit));
                    count2init(runinit) = 1;
                end
            elseif classresinit(runinit) == 2
                if udiffinit(runinit) > 0
                    count1init(runinit) = abs(udiffinit(runinit));
                    count2init(runinit) = 1;
                else
                    count1init(runinit) = 0;
                    count2init(runinit) = 0;
                end
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
            % find the index of observation with maximal mistake
            max_count1init = max(count1init);
            idx3init = find(count1init == max_count1init);
            
            % change this index
            if classresinit(idx3init) == 1
                classresinit(idx3init) = 2;
            else
                classresinit(idx3init) = 1;
            end
            
            % calcucalte the difference between U(with) and U(without) again
            meplus = mean(iinit(classresinit == 2));
            if isnan(meplus)
                meplus = 0;
            end
            meminus = mean(iinit(classresinit == 1));
            if isnan(meminus)
                meminus = 0;
            end
            mbplus = mean(dinit(classresinit == 2));
            if isnan(mbplus)
                mbplus = 0;
            end
            mbminus = mean(dinit(classresinit == 1));
            if isnan(mbminus)
                mbminus = 0;
            end
            veminus = std(iinit(classresinit == 1))*std(iinit(classresinit == 1));
            if isnan(veminus)
                veminus = 0;
            end
            vbplus = std(dinit(classresinit == 2))*std(dinit(classresinit == 2));
            if isnan(vbplus)
                vbplus = 0;
            end
            vbminus = std(dinit(classresinit == 1))*std(dinit(classresinit == 1));
            if isnan(vbminus)
                vbminus = 0;
            end
            pr = mean(classresinit)-1;
            proizvplus = mean(dinit(classresinit == 2))*mean(iinit(classresinit == 2));
            if isnan(proizvplus)
                proizvplus = 0;
            end
            proizvminus = mean(dinit(classresinit == 1))*mean(iinit(classresinit == 1));
            if isnan(proizvminus)
                proizvminus = 0;
            end
            % calculate mbplustwo
            mbplustwo = zeros(n_runsinit, 1);
            % iterate through each row in the table
            for idop = 1:n_runsinit
                % get the value from column 'i' in the current row
                valuei = iinit(idop);
                % get the stored indices for the current 'i' value
                indices = selectedIndices{idop};
                % get the values in column 'class' and 'd' for the selected indices
                selectedclass = classresinit(indices);
                selectedd = dinit(indices);
                % calculate the mean of values in column 'd' where 'class' == 2
                mbplustwodop = mean(selectedd(selectedclass == 2));
                if isnan(mbplustwodop)
                    mbplustwodop = 0;
                end
                % assign the mean to column 'mbplustwo' in the current row
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
            for runinit = 1:n_runsinit
                if classresinit(runinit) == 1
                    if udiffinit(runinit) >= 0
                        count1init(runinit) = 0;
                        count2init(runinit) = 0;
                    else
                        count1init(runinit) = abs(udiffinit(runinit));
                        count2init(runinit) = 1;
                    end
                elseif classresinit(runinit) == 2
                    if udiffinit(runinit) > 0
                        count1init(runinit) = abs(udiffinit(runinit));
                        count2init(runinit) = 1;
                    else
                        count1init(runinit) = 0;
                        count2init(runinit) = 0;
                    end
                end
            end
            
            % counter for count2init==1 - we will further escape the loop if we can't
            % progress to 0
            % and counter for the number of runs
            if (sum(count2init)==1)
                count_onesinit = count_onesinit+1;
            else
                count_allrunsinit = count_allrunsinit+1;
                if (count_allruns_storeinit(outerruninit)<count_allrunsinit)
                    count_allruns_storeinit(outerruninit) = count_allrunsinit;
                end
            end
            
            min_mistakeinit = min(sum(count2init), sum2_oldinit);
            
            % stop if no mistakes, or we loop on 1 mistake, or we hit the limit of 10000 runs
            % (in this case we are within the loop of the length n, say 17 12 5 3 9 17 ...
            % mistakes, so we stop on 3, also checking that this number of mistakes doesn't exceed
            % sqrt(n_runsinit)), or we hit the limit of 11000 runs (happens when
            % in a previous case we have a loop of say 4 4 4 4 4... mistakes)
            % practically for the smaller model count2init=0 all the
            % time, this block of conditions is needed only for the large
            % model (n_runs=10000)
            if (sum(count2init)<1)
                break;
            elseif (count_onesinit == 10)
                break;
            elseif (count_allrunsinit > 9999 && sum(count2init)>sum2_oldinit && min_mistakeinit<sqrt(n_runsinit))
                break;
            elseif (count_allrunsinit == 11000)
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
    
    % get the sorted indices in descending order
    [sortedValues, sortedIndices] = sort(sumudiffinit_values, 'ascend');
    
    % extract the indices of the half of the largest values
    indicesOfSmallest = sortedIndices(1:numSmallest);
    
    % copy the matrix
    reshapedMatrix1 = reshapedMatrix;
    
    % take half best runs of the small model and calculate the class for
    % each dot
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
    % observations from the small model are in the center of the grid
    % cells; then assign the class from the small model to all the dots in
    % the grid cell
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
    
    % store the generated values
    for run = 1:n_runs
        a_values(run) = a;
        e_values(run) = e;
        f_values(run) = f;
        g_values(run) = g;
        h_values(run) = h;
        j_values(run) = j;
        k_values(run) = k;
        l_values(run) = l;
        m_values(run) = m;
        n_values(run) = n;
        o_values(run) = o;
        odop_values(run) = odop;
    end
    
    % calcucalte the difference between U(with) and U(without)
    count_ones = 0;
    count_allruns = 0;
    meplus = mean(i(classres == 2));
    meminus = mean(i(classres == 1));
    mbplus = mean(d(classres == 2));
    mbminus = mean(d(classres == 1));
    veminus = std(i(classres == 1))*std(i(classres == 1));
    vbplus = std(d(classres == 2))*std(d(classres == 2));
    vbminus = std(d(classres == 1))*std(d(classres == 1));
    pr = mean(classres)-1;
    % calculate mbplustwo
    mbplustwo = zeros(n_runs, 1);
    selectedIndices = cell(n_runs, 1);
    % iterate through each row in the table
    for idop = 1:n_runs
        % get the value from column in the current row
        valuei = i(idop);
        % find the indices of the X closest values to the current value in column
        [~, indices] = mink(abs(i - valuei), 200);
        selectedIndices{idop} = indices;
        selectedclass = classres(indices);
        selectedd = d(indices);
        % count the number of 2s in column for the selected indices
        mbplustwodop = mean(selectedd(selectedclass == 2));
        if isnan(mbplustwodop)
            mbplustwodop = 0;
        end
        % assign the count to column in the current row
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
    for run = 1:n_runs
        if classres(run) == 1
            if udiff(run) >= 0
                count1(run) = 0;
                count2(run) = 0;
            else
                count1(run) = abs(udiff(run));
                count2(run) = 1;
            end
        elseif classres(run) == 2
            if udiff(run) > 0
                count1(run) = abs(udiff(run));
                count2(run) = 1;
            else
                count1(run) = 0;
                count2(run) = 0;
            end
        end
    end
    
    % store classres/count1/count2 (further will be used as the values on
    % the previous run)
    class_start = classres;
    count1_start = count1;
    count2_start = count2;
    
    sum1_old = sum(count1);
    sum2_old = sum(count2);
    
    % loop to reassign class to get rid of mistakingly assigned classes
    while true
        % find the index of observation with maximal mistake
        max_count1 = max(count1);
        idx3 = find(count1 == max_count1);
        
        % change this index
        if classres(idx3) == 1
            classres(idx3) = 2;
        else
            classres(idx3) = 1;
        end
        
        % calcucalte U(with) - U(without) with new classes
        meplus = mean(i(classres == 2));
        if isnan(meplus)
            meplus = 0;
        end
        meminus = mean(i(classres == 1));
        if isnan(meminus)
            meminus = 0;
        end
        mbplus = mean(d(classres == 2));
        if isnan(mbplus)
            mbplus = 0;
        end
        mbminus = mean(d(classres == 1));
        if isnan(mbminus)
            mbminus = 0;
        end
        veminus = std(i(classres == 1))*std(i(classres == 1));
        if isnan(veminus)
            veminus = 0;
        end
        vbplus = std(d(classres == 2))*std(d(classres == 2));
        if isnan(vbplus)
            vbplus = 0;
        end
        vbminus = std(d(classres == 1))*std(d(classres == 1));
        if isnan(vbminus)
            vbminus = 0;
        end
        pr = mean(classres)-1;
        proizvplus = mean(d(classres == 2))*mean(i(classres == 2));
        if isnan(proizvplus)
            proizvplus = 0;
        end
        proizvminus = mean(d(classres == 1))*mean(i(classres == 1));
        if isnan(proizvminus)
            proizvminus = 0;
        end
        % calculate mbplustwo
        mbplustwo = zeros(n_runs, 1);
        % iterate through each row in the table
        for idop = 1:n_runs
            % get the value from column 'i' in the current row
            valuei = i(idop);
            % get the stored indices for the current 'i' value
            indices = selectedIndices{idop};
            % get the values in column 'class' and 'd' for the selected indices
            selectedclass = classres(indices);
            selectedd = d(indices);
            % calculate the mean of values in column 'd' where 'class' == 2
            mbplustwodop = mean(selectedd(selectedclass == 2));
            if isnan(mbplustwodop)
                mbplustwodop = 0;
            end
            % assign the mean to column 'mbplustwo' in the current row
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
        for run = 1:n_runs
            if classres(run) == 1
                if udiff(run) >= 0
                    count1(run) = 0;
                    count2(run) = 0;
                else
                    count1(run) = abs(udiff(run));
                    count2(run) = 1;
                end
            elseif classres(run) == 2
                if udiff(run) > 0
                    count1(run) = abs(udiff(run));
                    count2(run) = 1;
                else
                    count1(run) = 0;
                    count2(run) = 0;
                end
            end
        end
        
        % counter for count2init==1 - we will further escape the loop if we can't
        % progress to 0
        % and counter for the number of runs
        if (sum(count2)==1)
            count_ones = count_ones+1;
        else
            count_allruns = count_allruns+1;
            if (count_allruns_store(outerrun)<count_allruns)
                count_allruns_store(outerrun) = count_allruns;
            end
        end
        
        min_mistake = min(sum(count2), sum2_old);
        
        % stop if no mistakes, or we looped on 1 mistake, or we hit the limit of 10000 runs
        % (in this case we are within the loop, say 17 12 5 3 9 17 ...
        % mistakes, so we stop on 3, also checking that this number of mistakes doesn't exceed
        % 30), or we hit the limit of 11000 runs (happens when
        % in a previous case we have a loop of say 4 4 4 4 4... mistakes)
        if (sum(count2)<1)
            break;
        elseif (count_ones == 10)
            break;
        elseif (count_allruns > 9999 && sum(count2)>sum2_old && min_mistake<30)
            break;
        elseif (count_allruns == 11000)
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
    
    % store draws in a table    
    second_table = [second_table; table(i,d,classres, udiff, count1, count2, sum1ch, sum2ch, ...
        class_start, count1_start, count2_start, a_values, e_values, f_values, g_values, h_values, ...
        j_values, k_values, l_values, m_values, n_values, o_values, odop_values, mbplustwo)];
    
    % we decide which classes to take - from the last run or from the second
    % to last, and then we store the results
    if (sum(count2)>sum2_old)
        store_eminus_check = mean(second_table.i(second_table.class_start==1));
        store_eplus_check = mean(second_table.i(second_table.class_start==2));
        store_bminus_check = mean(second_table.d(second_table.class_start==1));
        store_bplus_check = mean(second_table.d(second_table.class_start==2));
        store_pr_check = mean(second_table.class_start)-1;
        store_veminus_check = std(second_table.i(second_table.class_start==1))*std(second_table.i(second_table.class_start==1));
        store_vbminus_check = std(second_table.d(second_table.class_start==1))*std(second_table.d(second_table.class_start==1));
        store_vbplus_check = std(second_table.d(second_table.class_start==2))*std(second_table.d(second_table.class_start==2));
        store_proizvminus_check = mean(second_table.d(second_table.class_start==1))*mean(second_table.i(second_table.class_start==1));
        store_proizvplus_check = mean(second_table.d(second_table.class_start==2))*mean(second_table.i(second_table.class_start==2)); 
        store_bplustwo_dop2_check = zeros(10000, 1);
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
            store_bplustwo_dop_check = mean(second_table.d(second_table.class_start==2 & ...
                second_table.i>xedges3(min_bound) & second_table.i<xedges3(max_bound)));
            if isnan(store_bplustwo_dop_check)
                store_bplustwo_dop_check = 0;
            end
            store_bplustwo_dop2_check(rungraph7) = store_bplustwo_dop_check;
        end
        store_bplustwo_check = smoothdata(store_bplustwo_dop2_check, 'gaussian', 1000);
        store_rp1_check = a*(n^2)*(store_pr_check*store_bplus_check*m...
            +(1-store_pr_check)*(store_veminus_check+m*store_bminus_check))...
            +a*(n^2)*store_pr_check*(1-store_pr_check)*((store_eplus_check-store_eminus_check)^2)...
            +0.5*(a^3)*(n^4)*store_pr_check*(1-store_pr_check)*((store_bplus_check*m-store_veminus_check-store_bminus_check*m)^2)...
            -1.5*(a^2)*(n^3)*store_pr_check*m*store_proizvplus_check...
            -1.5*(a^2)*(n^3)*(1-store_pr_check)*(store_eminus_check*store_veminus_check+m*store_proizvminus_check)...
            +1.5*(a^2)*(n^3)*(store_pr_check*store_eplus_check...
            +(1-store_pr_check)*store_eminus_check)*(store_pr_check*store_bplus_check*m...
            +(1-store_pr_check)*(store_veminus_check+store_bminus_check*m));
    else
        store_eminus_check = mean(second_table.i(second_table.classres==1));
        store_eplus_check = mean(second_table.i(second_table.classres==2));
        store_bminus_check = mean(second_table.d(second_table.classres==1));
        store_bplus_check = mean(second_table.d(second_table.classres==2));
        store_pr_check = mean(second_table.classres)-1;
        store_veminus_check = std(second_table.i(second_table.classres==1))*std(second_table.i(second_table.classres==1));
        store_vbminus_check = std(second_table.d(second_table.classres==1))*std(second_table.d(second_table.classres==1));
        store_vbplus_check = std(second_table.d(second_table.classres==2))*std(second_table.d(second_table.classres==2));
        store_proizvminus_check = mean(second_table.d(second_table.classres==1))*mean(second_table.i(second_table.classres==1));
        store_proizvplus_check = mean(second_table.d(second_table.classres==2))*mean(second_table.i(second_table.classres==2));
        store_bplustwo_dop2_check = zeros(10000, 1);
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
            store_bplustwo_dop_check = mean(second_table.d(second_table.classres==2 & ...
                second_table.i>xedges3(min_bound) & second_table.i<xedges3(max_bound)));
            if isnan(store_bplustwo_dop_check)
                store_bplustwo_dop_check = 0;
            end
            store_bplustwo_dop2_check(rungraph7) = store_bplustwo_dop_check;
        end
        store_bplustwo_check = smoothdata(store_bplustwo_dop2_check, 'gaussian', 1000);
        store_rp1_check = a*(n^2)*(store_pr_check*store_bplus_check*m+(1-store_pr_check)*(store_veminus_check+m*store_bminus_check))...
            +a*(n^2)*store_pr_check*(1-store_pr_check)*((store_eplus_check-store_eminus_check)^2)...
            +0.5*(a^3)*(n^4)*store_pr_check*(1-store_pr_check)*((store_bplus_check*m-store_veminus_check-store_bminus_check*m)^2)...
            -1.5*(a^2)*(n^3)*store_pr_check*m*store_proizvplus_check...
            -1.5*(a^2)*(n^3)*(1-store_pr_check)*(store_eminus_check*store_veminus_check+m*store_proizvminus_check)...
            +1.5*(a^2)*(n^3)*(store_pr_check*store_eplus_check...
            +(1-store_pr_check)*store_eminus_check)*(store_pr_check*store_bplus_check*m...
            +(1-store_pr_check)*(store_veminus_check+store_bminus_check*m));
    end
    
    store_eminus = store_eminus_check;
    store_eplus = store_eplus_check;
    store_bminus = store_bminus_check;
    store_bplus = store_bplus_check;
    store_pr = store_pr_check;
    store_veminus = store_veminus_check;
    store_vbminus = store_vbminus_check;
    store_vbplus = store_vbplus_check;
    store_proizvminus = store_proizvminus_check;
    store_proizvplus = store_proizvplus_check;
    store_bplustwo = store_bplustwo_check;
    store_rp1 = store_rp1_check;
    
    % set the 10000*10000 grid for the second step of the algorithm
    dgrid2 = linspace(0,2,10001);
    igrid2 = icdf(p_i, linspace(0,1,10001));
    [X2,Y2] = meshgrid(igrid2,dgrid2);
    xedges2 = igrid2;
    % replace infinity
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
    
    % calcucalte U(with) - U(without) for each dot and store 
    % 0 or 1 for the class and Udiff itself
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
    % calculate mean and variance of Udiff
    mean_final_grid2 = mean(final_grid2(:));
    var_final_grid2 = var(final_grid2(:));
    % count the number of ones in the matrix (when we intervene)
    num_ones = sum(final_grid(:) == 1);
    % store if we've faced a loop
    if (count_allruns_store(outerrun) > 9999)
        infin_store=1;
    else
        infin_store=0;
    end
    % calculate the proportion of ones - Prop^+
    proportion_ones = num_ones / numel(final_grid);
    % calculate the number of zeros and ones for some \beta
    num_zeros10 = sum(final_grid(:,10000) == 0);
    num_zeros9 = sum(final_grid(:,9000) == 0);
    num_zeros8 = sum(final_grid(:,8000) == 0);
    num_zeros7 = sum(final_grid(:,7000) == 0);
    num_zeros6 = sum(final_grid(:,6000) == 0);
    num_zeros5 = sum(final_grid(:,5000) == 0);
    num_zeros4 = sum(final_grid(:,4000) == 0);
    num_zeros3 = sum(final_grid(:,3000) == 0);
    num_zeros2 = sum(final_grid(:,2000) == 0);
    num_zeros1 = sum(final_grid(:,1000) == 0);
    num_zeros0 = sum(final_grid(:,1) == 0);
    num_ones10 = sum(final_grid(:,10000) == 1);
    num_ones9 = sum(final_grid(:,9000) == 1);
    num_ones8 = sum(final_grid(:,8000) == 1);
    num_ones7 = sum(final_grid(:,7000) == 1);
    num_ones6 = sum(final_grid(:,6000) == 1);
    num_ones5 = sum(final_grid(:,5000) == 1);
    num_ones4 = sum(final_grid(:,4000) == 1);
    num_ones3 = sum(final_grid(:,3000) == 1);
    num_ones2 = sum(final_grid(:,2000) == 1);
    num_ones1 = sum(final_grid(:,1000) == 1);
    num_ones0 = sum(final_grid(:,1) == 1);
    % calculate the number of zeros and ones for some \varepsilon_2
    onesinminustwosigma = sum(final_grid(228,:) == 1);
    propminustwosigma = onesinminustwosigma / size(final_grid,2);
    onesinminussigma = sum(final_grid(1587,:) == 1);
    propminussigma = onesinminussigma / size(final_grid,2);
    onesinminushalfsigma = sum(final_grid(3086,:) == 1);
    propminushalfsigma = onesinminushalfsigma / size(final_grid,2);
    onesinminusquartersigma = sum(final_grid(4013,:) == 1);
    propminusquartersigma = onesinminusquartersigma / size(final_grid,2);
    onesinminuszerosigma = sum(final_grid(5000,:) == 1);
    onesinminuszerosigma1 = sum(final_grid(5001,:) == 1);
    propminuszerosigma = (onesinminuszerosigma+onesinminuszerosigma1)/(2*size(final_grid,2));
    onesinplusquartersigma = sum(final_grid(5988,:) == 1);
    propplusquartersigma = onesinplusquartersigma / size(final_grid,2);
    onesinplushalfsigma = sum(final_grid(6915,:) == 1);
    propplushalfsigma = onesinplushalfsigma / size(final_grid,2);
    onesinplussigma = sum(final_grid(8414,:) == 1);
    propplussigma = onesinplussigma / size(final_grid,2);
    onesinplustwosigma = sum(final_grid(9773,:) == 1);
    propplustwosigma = onesinplustwosigma / size(final_grid,2);
    percentiles_eminus = 10000*normcdf(store_eminus, mui, sigmai);
    percentiles_eplus = 10000*normcdf(store_eplus, mui, sigmai);

    % store all the equilibrium characteristics for the case Prop^+==1
    if (num_ones == 100000000)
        % find the indices of the zero/one elements in the matrix
        zero_indices = NaN;
        one_indices = find(final_grid == 1);
        % find the column index of the first and the last zero/one element
        min_i = NaN;
        max_i = NaN;
        min_d = NaN;
        max_d = NaN;
        min_i1 = min(mod(one_indices - 1, size(final_grid, 1)) + 1);
        max_i1 = max(mod(one_indices - 1, size(final_grid, 1)) + 1);
        min_d1 = min(floor((one_indices - 1) / size(final_grid, 1)) + 1);
        max_d1 = max(floor((one_indices - 1) / size(final_grid, 1)) + 1);
        % compute the row indices of the zero elements
        row_indices = NaN;
        % calculate the average row index
        avg_row = NaN;
        % compute the column indices of the zero elements
        col_indices = NaN;
        % calculate the average column index
        avg_col = NaN;
        % compute the row indices of the one elements
        row_indices1 = floor((one_indices - 1) / size(final_grid, 1)) + 1;
        % calculate the average row index
        avg_row1 = mean(row_indices1);
        % compute the column indices of the one elements
        col_indices1 = mod(one_indices - 1, size(final_grid, 1)) + 1;
        % calculate the average column index
        avg_col1 = mean(col_indices1);
        % calculate the variance of the indices
        variance_d_indices = NaN;
        variance_d_indices1 = var(row_indices1);
        variance_i_indices = NaN;
        variance_i_indices1 = var(col_indices1);
        % add another line to the table with all the draws
        table4 = [table4; table(proportion_ones, min_i, max_i, min_d, max_d, min_i1, max_i1, ...
            min_d1, max_d1, variance_d_indices, variance_d_indices1, variance_i_indices, ...
            variance_i_indices1, num_zeros10, num_zeros9, num_zeros8, num_zeros7, num_zeros6, ...
            num_zeros5, num_zeros4, num_zeros3, num_zeros2, num_zeros1, num_zeros0, num_ones10, ...
            num_ones9, num_ones8, num_ones7, num_ones6, num_ones5, num_ones4, num_ones3, ...
            num_ones2, num_ones1, num_ones0, avg_row, avg_col, avg_row1, avg_col1, ...
            mean_final_grid2, var_final_grid2, propminustwosigma, propminussigma, ...
            propminushalfsigma, propminusquartersigma, propminuszerosigma, propplusquartersigma, ...
            propplushalfsigma, propplussigma, propplustwosigma, a, e, f, g, h, j, k, l, m, n, o, ...
            odop, mui, sigmai, store_eminus, store_eplus, store_bminus, store_bplus, store_pr, ...
            store_veminus, store_vbminus, store_vbplus, store_proizvminus, store_proizvplus, ...
            store_rp1, percentiles_eminus, percentiles_eplus, infin_store, min_mistake, count_allruns)];
        % store all the equilibrium characteristics for the case Prop^+==0
    elseif (num_ones == 0)
        % find the indices of the zero/one elements in the matrix
        zero_indices = find(final_grid == 0);
        one_indices = NaN;
        % find the column index of the first and the last zero/one element
        min_i = min(mod(zero_indices - 1, size(final_grid, 1)) + 1);
        max_i = max(mod(zero_indices - 1, size(final_grid, 1)) + 1);
        min_d = min(floor((zero_indices - 1) / size(final_grid, 1)) + 1);
        max_d = max(floor((zero_indices - 1) / size(final_grid, 1)) + 1);
        min_i1 = NaN;
        max_i1 = NaN;
        min_d1 = NaN;
        max_d1 = NaN;
        % compute the row indices of the zero elements
        row_indices = floor((zero_indices - 1) / size(final_grid, 1)) + 1;
        % calculate the average row index
        avg_row = mean(row_indices);
        % compute the column indices of the zero elements
        col_indices = mod(zero_indices - 1, size(final_grid, 1)) + 1;
        % calculate the average column index
        avg_col = mean(col_indices);
        % compute the row indices of the one elements
        row_indices1 = NaN;
        % calculate the average row index
        avg_row1 = NaN;
        % compute the column indices of the one elements
        col_indices1 = NaN;
        % calculate the average column index
        avg_col1 = NaN;
        % calculate the variance of the indices
        variance_d_indices = var(row_indices);
        variance_d_indices1 = NaN;
        variance_i_indices = var(col_indices);
        variance_i_indices1 = NaN;
        % add another line to the table with all the draws
        table4 = [table4; table(proportion_ones, min_i, max_i, min_d, max_d, min_i1, max_i1, ...
            min_d1, max_d1, variance_d_indices, variance_d_indices1, variance_i_indices, ...
            variance_i_indices1, num_zeros10, num_zeros9, num_zeros8, num_zeros7, num_zeros6, ...
            num_zeros5, num_zeros4, num_zeros3, num_zeros2, num_zeros1, num_zeros0, num_ones10, ...
            num_ones9, num_ones8, num_ones7, num_ones6, num_ones5, num_ones4, num_ones3, ...
            num_ones2, num_ones1, num_ones0, avg_row, avg_col, avg_row1, avg_col1, ...
            mean_final_grid2, var_final_grid2, propminustwosigma, propminussigma, ...
            propminushalfsigma, propminusquartersigma, propminuszerosigma, ...
            propplusquartersigma, propplushalfsigma, propplussigma, propplustwosigma, ...
            a, e, f, g, h, j, k, l, m, n, o, odop, mui, sigmai, store_eminus, ...
            store_eplus, store_bminus, store_bplus, store_pr, store_veminus, ...
            store_vbminus, store_vbplus, store_proizvminus, store_proizvplus, ...
            store_rp1, percentiles_eminus, percentiles_eplus, infin_store, min_mistake, count_allruns)];
        % store all the equilibrium characteristics for Prop^+ != 0;1
    else
        % find the indices of the zero/one elements in the matrix
        zero_indices = find(final_grid == 0);
        one_indices = find(final_grid == 1);
        % find the column index of the first and the last zero/one element
        min_i = min(mod(zero_indices - 1, size(final_grid, 1)) + 1);
        max_i = max(mod(zero_indices - 1, size(final_grid, 1)) + 1);
        min_d = min(floor((zero_indices - 1) / size(final_grid, 1)) + 1);
        max_d = max(floor((zero_indices - 1) / size(final_grid, 1)) + 1);
        min_i1 = min(mod(one_indices - 1, size(final_grid, 1)) + 1);
        max_i1 = max(mod(one_indices - 1, size(final_grid, 1)) + 1);
        min_d1 = min(floor((one_indices - 1) / size(final_grid, 1)) + 1);
        max_d1 = max(floor((one_indices - 1) / size(final_grid, 1)) + 1);
        % compute the row indices of the zero elements
        row_indices = floor((zero_indices - 1) / size(final_grid, 1)) + 1;
        % calculate the average row index
        avg_row = mean(row_indices);
        % compute the column indices of the zero elements
        col_indices = mod(zero_indices - 1, size(final_grid, 1)) + 1;
        % calculate the average column index
        avg_col = mean(col_indices);
        % compute the row indices of the one elements
        row_indices1 = floor((one_indices - 1) / size(final_grid, 1)) + 1;
        % calculate the average row index
        avg_row1 = mean(row_indices1);
        % compute the column indices of the one elements
        col_indices1 = mod(one_indices - 1, size(final_grid, 1)) + 1;
        % calculate the average column index
        avg_col1 = mean(col_indices1);
        % calculate the variance of the indices
        variance_d_indices = var(row_indices);
        variance_d_indices1 = var(row_indices1);
        variance_i_indices = var(col_indices);
        variance_i_indices1 = var(col_indices1);
        % add another line to the table with all the draws
        table4 = [table4; table(proportion_ones, min_i, max_i, min_d, max_d, min_i1, ...
            max_i1, min_d1, max_d1, variance_d_indices, variance_d_indices1, ...
            variance_i_indices, variance_i_indices1, num_zeros10, num_zeros9, ...
            num_zeros8, num_zeros7, num_zeros6, num_zeros5, num_zeros4, num_zeros3, ...
            num_zeros2, num_zeros1, num_zeros0, num_ones10, num_ones9, num_ones8, ...
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