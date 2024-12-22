
%% Houskeeping

% HS_vals - Short for Hot-Spots values (or scores) matrix. 
%           These matrices pack together hot spot valus for a given group of dimers (or chains) where each column corresponds to a spasific dimer (or chain), 
%                                                                                          e.g. HS_vals(:,4) are the ascores for all hot spots calcualte for dimer No. 4
%           In most of the the code variabes of HS_vals tyoe are packed into cell arrays,  e.g., HS_vals_Chains_cell) 
%
%           size(hs_vasl) should be: [number of hot spots, number of dimers (or chains)] => in most casess [~18,12] or [~18,24]          

[z, z, maps_scores_list] = xlsread('maps_scores_list_06_Nov_23.csv'); 

maps_scores_list  = maps_scores_list(2:end,:)    ;
maps_names=maps_scores_list(:,4);
num_maps_types_from_list=size(maps_scores_list,1);
num_maps_types_tot      =num_maps_types_from_list; % This var is defined to allow adding values of symulated maps (in thoes case it will be increased accordingly)

[z, z, dublets_list] = xlsread('dublets_list.csv'); dublets_list= dublets_list(2:end,:); % A list of chain pairs (monomers in dimers) each line documets the two associated chains
dimers_couplings = dublets_list(:,[3,6]);
dimers_names     = dublets_list(:,10);
chains_inds      = [[dublets_list{:,7}]',[dublets_list{:,8}]'];
num_dimers=size(dublets_list,1);

[z, z, QS_res_list]=xlsread('Q_Scores_Res_Coupling_FN.csv'); QS_res_list= QS_res_list(2:end,:); % A list of residus numbers of Bfr1 and Bfr2 that their Q_scores are compared


res_Bfr1=[QS_res_list{:,1}]';
res_Bfr2=[QS_res_list{:,2}]';
res_include=[QS_res_list{:,3}]'; 
res_include=res_include>0;
num_res_Q_scores=sum(res_include>0);

base_dir='C:\Users\BFR12\';


%% Building: Q-scores matrices
%
% This segment will iterate on every normlization method of the data and on every map type (i.e. sorecers for Q-scores) 
%      and coalate the tables of Q-scores for each dimer or chain per map and normalization method
%
maps_names=maps_scores_list(:,1);
num_norm_types=1; % Number of diffretn ways for pre-processing the Q-scores 
% Each matrix in QS_vals_dimers_Cell corresponsd to the hot spots values (or scores) calcualted for the correspoding map in maps_cell (named in QS_notes)
%
% The column of the matrix correspond to the dimers i.e. the column vector QS_vals_dimers_Cell{3,1}(:,1} holds the Q-scores of the dimer AQ (the must be exactly 12 columns)
QS_vals_dimers_Cell = repmat({nan([num_res_Q_scores, num_dimers]  )}, [num_norm_types, num_maps_types_tot]); %
QS_vals_chains_Cell = repmat({nan([num_res_Q_scores, num_dimers*2])}, [num_norm_types, num_maps_types_tot]); %
QS_notes            = repmat({'empty-empty-empty'}             , [num_norm_types, num_maps_types_tot]);

res_num_list=[];
for q1=1:num_maps_types_tot
    [z, z, Q_scores_Bfr1_q1]=xlsread([base_dir, maps_scores_list{q1,2}]); Q_scores_Bfr1_q1 = Q_scores_Bfr1_q1(2:end,:);
    [z, z, Q_scores_Bfr2_q1]=xlsread([base_dir, maps_scores_list{q1,3}]); Q_scores_Bfr2_q1 = Q_scores_Bfr2_q1(2:end,:);
    
    for q2=1:num_dimers
        chain_1_q2=dimers_couplings{q2,1};
        chain_2_q2=dimers_couplings{q2,2};
        
        for q3=1:num_res_Q_scores %
            res_Bfr1_q3=res_Bfr1(q3);
            res_Bfr2_q3=res_Bfr2(q3);
            res_include_q3=res_include(q3);

            [Q_score_val_Bfr1_q123_c1, ind]=find_Q_score_val(Q_scores_Bfr1_q1, chain_1_q2, res_Bfr1_q3*res_include_q3);
            [Q_score_val_Bfr1_q123_c2, ind]=find_Q_score_val(Q_scores_Bfr1_q1, chain_2_q2, res_Bfr1_q3*res_include_q3);
            [Q_score_val_Bfr2_q123_c1, ind]=find_Q_score_val(Q_scores_Bfr2_q1, chain_1_q2, res_Bfr2_q3*res_include_q3);
            [Q_score_val_Bfr2_q123_c2, ind]=find_Q_score_val(Q_scores_Bfr2_q1, chain_2_q2, res_Bfr2_q3*res_include_q3);   

            % Iterating over q1, q2, q3 
                           % q1 - maps
                           % q2 - dimers and chains
                           % q3 - res numbers of Hot-Spots 
            norm_type=1; 
            QS_vals_dimers_Cell{norm_type, q1}(q3,q2) = (Q_score_val_Bfr1_q123_c1+Q_score_val_Bfr1_q123_c2) - (Q_score_val_Bfr2_q123_c1+Q_score_val_Bfr2_q123_c2);
            % NOTE: The current way th eprogram is run the '-' sign is irrelvant as the terms befor and after the sign are never non zero at toghter
            if q3==1; QS_notes{norm_type, q1}= [maps_names{q1},' delta of sums']; end

        end % q3 - residue
    end
end % q1 - maps types






% %% Mean analysis
% %
% % If the different Hot spot values are roughly on the same scale => that the mean of the HS-vals of each dimer (chain) can teach us a lot on its classfication 
% %
% % Iterating over the Hot spot values of dimers calcualted for all 'Hot-Spots-Mask-Norm' maps { i.e., iterating over maps_cell(3,q1) } 
% % and calculateing the column-wise mean.
% %  The result of the calculation is a 12 elements row vector, where each element is the mean of all the hot spots values of a spasific dimer
% 
% map_norm_type=1;
% for q1=1:num_maps_types_tot
%     mean_QS_q1    =mean(QS_vals_dimers_Cell{map_norm_type,q1},1);
%     mean_mat(q1,:)=mean_QS_q1;
%     if q1==1 
%         [h,ind_srt,srt_names]=plot_sorted_items(mean_QS_q1', dimers_names, ['Mean of hot spot scores of all dimers']);
%         %set(h,'name',['TH sig=',num2str(sig)]);
% 
%         hold on
%         legend(QS_notes{map_norm_type,q1});
%     else
% 
%         plot([1:numel(mean_QS_q1)]', mean_QS_q1(ind_srt),'--o' , "DisplayName", QS_notes{map_norm_type,q1});
%         %text([1:numel(mean_HS_q1)]', mean_HS_q1(ind_srt), srt_names, 'VerticalAlignment', 'top'); 
%         %legend(QS_notes{3,q1});
%     end
% end
% hold off 

%% PCA
map_norm_type=1;
X1=[1:12]';
for q1=1:num_maps_types_tot
    % Centreing the Hot-Spots values
    % i.e. substracting the mean of each Hot-Spot (mean clacualted along the rows) from all vectors
    % 
    
    if q1==1 % Calculating the singular vectors of the refrence group {in this case the HS_vals of the Sym_relaxed map after HotSpot mask normalization} 
        centering_vect = mean(QS_vals_dimers_Cell{map_norm_type,q1},2); % Estimating the shift from the center of each hot-spot scors 
        centerd_HS_vals_q1=QS_vals_dimers_Cell{map_norm_type,q1} - centering_vect;
        
        
        [coeff,score,~,~,explained,~]=pca_of_HS_vals(centerd_HS_vals_q1); % To understad the directioanltey of the feature vectors read the doc of the function pca_of_HS_vals()
        
        % Calculating the projection of each hot-spots values vector on the 1st principal component vector i.e. the projection of the columns of HS_vals on coeff(:,1)
        % Note that these values are the same as in score(:,1)
        % But can be also calcualted by: 
        % proj_on_1st_PC_q1=(coeff(:,1)'*centerd_HS_vals_q1)';  <=== This will become important later on
        proj_on_1st_PC_q1=score(:,1); % By deff - the 1st column of score is the same as above!
        [h,ind_srt,srt_names]=plot_sorted_items(proj_on_1st_PC_q1, dimers_names, ['proj on 1st Principal Component (explained variance)=', num2str(explained(1))]);
        %set(h,'name',['TH sig=',num2str(sig)]);
        YMatrix1(:,q1)=proj_on_1st_PC_q1(ind_srt);
        hold on
        legend(QS_notes{map_norm_type,q1});   
    else 
        centerd_HS_vals_q1=QS_vals_dimers_Cell{map_norm_type,q1} - centering_vect; % Centering all HS_val vectors based on the centering vector calculated for the refrence group 
         % Calculating the projection of each hot-spots values vector on the 1st principal component vector i.e. the projection of the columns of HS_vals on coeff(:,1)
        proj_on_1st_PC_q1=(coeff(:,1)'*centerd_HS_vals_q1)';    
        plot([1:numel(proj_on_1st_PC_q1)]', proj_on_1st_PC_q1(ind_srt),'--o' , "DisplayName", QS_notes{map_norm_type,q1});
        YMatrix1(:,q1)=proj_on_1st_PC_q1(ind_srt);
    end
end
hold off 


%% SVD analysis 

map_norm_type=1;
for q1=1:num_maps_types_tot
    % Centreing the Hot-Spots values
    % i.e. substracting the mean of each Hot-Spot (mean clacualted along the rows) from all vectors

    if q1==1 % Calculating the singular vectors of the refrence group {in this case the HS_vals of the Sym_relaxed map after HotSpot mask normalization} 
        centering_vect    =mean(QS_vals_dimers_Cell{map_norm_type,q1},2); % Estimating the shift from the center of each hot-spot scors 
        centerd_HS_vals_q1=QS_vals_dimers_Cell{map_norm_type,q1} - centering_vect;
        [U,S,~]=svd(centerd_HS_vals_q1,"vector");  % Calcualting SVD (treating the columns of centerd_HS_vals_q1 (i.e. the list of values for each diemer) as the objects
                                                   % The  columns of U are the singualr vectors (U is an orthogonal matrix)
         U1 =U; U1(:,2:end)=0;                     % U1   keeps only 1st singular vector  
         U1z=U; U1z(:,1)   =0;                     % Uz1  keeps all the other singular vector 

        % *UC* Calculating the projection of each hot-spots values vector on the 1st singualr vector i.e. U(:,1) 
        proj_1st_SV_q1     =-U(:,1)'*centerd_HS_vals_q1         ;   % projection of all dimers on 1st Singualr Vector
        proj_2nd_SV_q1     =-U(:,2)'*centerd_HS_vals_q1         ;   % projection of all dimers on 2nd Singualr Vector
        proj_3rd_SV_q1     =-U(:,3)'*centerd_HS_vals_q1         ;   % projection of all dimers on 3rd Singualr Vector 
        Err_proj_1st_SV_q1 =vecnorm(U1z'*centerd_HS_vals_q1,2,1)   % norm of the residual HS_vectors after removing the projection above
        
        %sqrt_var_res_of_proj_q1=norm(proj_on_1st_simulated_SingualrVector_q1-centerd_HS_vals_q1);

        [h,ind_srt,srt_names]=plot_sorted_items_err_bar(proj_1st_SV_q1(:), Err_proj_1st_SV_q1/2, dimers_names, ['proj on 1st Singular Vector (explained variance)=', num2str(S(1).^2/sum(S.^2)*100)]);
        %set(h,'name',['TH sig=',num2str(sig)]);
 
        hold on
        legend(QS_notes{map_norm_type,q1}); 
        %plot(proj_2nd_SV_q1(ind_srt),'-s');

        figure % 2D plot SV1 and SV2
        plot(proj_1st_SV_q1(ind_srt),proj_2nd_SV_q1(ind_srt), '*')
        text(proj_1st_SV_q1(ind_srt),proj_2nd_SV_q1(ind_srt),dimers_names(ind_srt))
        xlabel('SV1');
        ylabel('SV2');
        title('SV1 vs SV2')

        figure % Scree plot
        plot(S/S(1), '-*')
        xlabel('Component Number');
        ylabel('Norm Eigen val');
        title('Scree plot')


        Z_rlx_vals=proj_1st_SV_q1(ind_srt)'
        Z_rlx_Err =Err_proj_1st_SV_q1(ind_srt)'
        Z_rlx_names_after_srt=dimers_names(ind_srt)


    else 
        centerd_HS_vals_q1=QS_vals_dimers_Cell{map_norm_type,q1} - centering_vect; % Centering all HS_val vectors based on the centering vector calculated for the refrence group 
        proj_1st_SV_q1    =-U(:,1)'*centerd_HS_vals_q1 ;          % projection of all dimers on 1st Singualr Vector 
        Err_proj_1st_SV_q1=vecnorm(U1z'*centerd_HS_vals_q1,2,1); % norm of the residual HS_vectors after removing the projection above
        % *UC* sqrt_var_res_of_proj_q1=norm(proj_on_1st_simulated_SingualrVector_q1-centerd_HS_vals_q1);
        figure(h)
        errorbar([1:numel(proj_1st_SV_q1)]', proj_1st_SV_q1(ind_srt), Err_proj_1st_SV_q1(ind_srt)/2, '--o' , "DisplayName", QS_notes{map_norm_type,q1});
    
        Z_sym_vals=proj_1st_SV_q1(ind_srt)'
        Z_sym_Err =Err_proj_1st_SV_q1(ind_srt)'
        

    
    end

end
hold off 





%% end sounds
beep
beep





%% FUNCS
function [Q_score_val, ind]=find_Q_score_val(Q_scores_list, chain_id, res_number)
% This function accept
% The searches the Q_scores_list for the relevent Q-score i.e. a score with chain id letter chani_id and residue number res_num and returen it and its index
% If no such Q_score is found it returens 0
% If more then one the function will rise and Error
%
% Call:
%[Q_score_val_Bfr1_q123_c2, ind]=find_Q_score_val(Q_scores_Bfr1_q1, chain_2_q2, res_Bfr1_q3);

            ind=find([Q_scores_list{1:end,1}]'==res_number & strcmpi(Q_scores_list(1:end,2),chain_id));
            switch  numel(ind)
                case 0
                    Q_score_val=0; % There is no Q score value for this chain
                case 1
                    Q_score_val=Q_scores_list{ind,3};
                otherwise
                    error('Error. \n The number of Q_Scores indices must be zero or one, got %d', ind)
            end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [blob_map_OUT, blob_size_OUT]=bigest_intersected_blob(logical_foot_print_1, logical_foot_print_2)
% This function accepts two same size logical 3D arrays of blobs (logical_foot_print_2 and logical_foot_print_2) which are logical footprints of 3D maps
% The function calcualte the intresection of the logical arrays (i.e. all places both logical foot pritns are true)
% The the function finds the largest connceted region in 3D (blob) in the intresected map
%
% Input
% logical_foot_print_1 - 3D logical 3D array (size [m,n,p]): to be intresected with the next input
% logical_foot_print_2 - 3D logical 3D array (size [m,n,p]): to be intresected with the previuse input
%
% Output
%
% blob_map_OUT  - 3D logical 3D array (size [m,n,p]): the largest intersected blob 
% blob_size_OUT - scalar : the number of true voxels in the map blob_map_OUT
%
% If the interasection between logical_foot_prints 1 & 2 is empty:
%                                                                  blob_map_OUT=false(size(logical_foot_print_1)) 
%                                                                  blob_size_OUT=0

blob_map_OUT=logical_foot_print_1 & logical_foot_print_2                              ;
blob_size_OUT=0                                                                       ;
        if any(blob_map_OUT(:))    
            blobs_intersects_struct=bwconncomp(blob_map_OUT,26)                       ;
            HS_numPixels_sect = cellfun(@numel,blobs_intersects_struct.PixelIdxList)  ;
            [blob_size_OUT,idx] = max(HS_numPixels_sect)                              ;
            blob_map_OUT=false(size(blob_map_OUT))                                    ;
            blob_map_OUT(blobs_intersects_struct.PixelIdxList{idx})=true              ;
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,ind_srt,srt_names]=plot_sorted_items(V_in, names_cell, ttl)
% sort the values in V_in
% applay the same permutatin on the cell array pf names names_cell
% plot the sorted values accoriding with x coordiante as index
% print the name near each value point on the graph

% o="#D95319"; % Orange [0.8500 0.3250 0.0980]
% b="#0072BD"; % Blue   [0 0.4470 0.7410]
% x="#7E2F8E"; % Purp   [0.4940 0.1840 0.5560]

% abc =['A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X'];
% clrs={ o ; b ; b ; o ; b ; x ; o ; o ; b ; b ; o ; b ; b ; o ; b ; b ; o ; x ; b ; b ; o ; b ; b ; b };
if nargin==2
    ttl='';
end

[srt_V,ind_srt]=sort(V_in(:));
%srt_abc=[repmat([' '], [24,1]),  abc(ind_V)];
srt_names=names_cell(ind_srt);

h=figure('name','TH' );
plot([1:numel(V_in)]',srt_V,'-*','LineWidth',1,'Color', 'b');
t=text([1:numel(V_in)]',srt_V,srt_names,'VerticalAlignment','top');

%[t.Color]=clrs{ind_V};
title(ttl);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,ind_srt,srt_names]=plot_sorted_items_err_bar(V_in, Err_in, names_cell, ttl)
% sort the values in V_in
% applay the same permutatin on the vector of error bars: Err_in, and the cell array of names: names_cell
% plot the sorted values accoriding with x coordiante as index
% print the name near each value point on the graph

if nargin==3
    ttl='';
end

[srt_V,ind_srt]=sort(V_in(:));
%srt_abc=[repmat([' '], [24,1]),  abc(ind_V)];

srt_Err  = Err_in(ind_srt);
srt_names= names_cell(ind_srt);

h=figure('name','TH');
errorbar([1:numel(V_in)]', srt_V, srt_Err, '-*', 'LineWidth',1,'Color', 'b');
t=text([1:numel(V_in)]',srt_V,srt_names,'VerticalAlignment','top');

title(ttl);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [coeff,score,latent,tsquared,explained,mu]=pca_of_HS_vals(Hot_spots_vals_in)
% See explantion below for the transpose operation of the input 
%
% Hot_spots_vals_in - [p,n] matrix of hot spots values calculated for each dimer or chain 
%                           each of the n columns corespond to a chain (or a dimer)
%                           each of the p rows are list of reuslts of one of the hot-spots calculated for eacho of the dimers (or chains)
%
% PCA function of matlab treats the rows of the matrix as the obseravtions, i.e., each row should correspond to the valuse measured for a spasific dimer (or chain)
% However, all HS_val matrices are ordered so each dimer (or chain) as columns NOT as roes thus they all HS_val matrices must  be "rotated" on their sides before being input into the PCA function (as stated above)

[coeff,score,latent,tsquared,explained,mu]=pca(Hot_spots_vals_in.'); % Note the .' operation!

% coeff = pca(X) returns the Principal Component Coefficients, also known as loadings, for the n-by-p data matrix X. 
% Input X [n,p] matrix:
%          the n Rows    of X correspond to a vector of observations 
%                        => i.e., group of features measured for a spesific object,
%                                 a realization of a multivariate experiment using several dectors,
%                                 the list of hot-spot scores of a spasific chain
%              p Columns of X correspond to variable measured in each observation 
%                        => i.e., list of values of a spesific feature measured from all objects, 
%                                 list of values from a spasifc detector,
%                                 list of scores of a spasific Hot-Spot)
% Output
%        coeff  The coefficient matrix is p-by-p. 
%               Each column of coeff is a "component" i.e., one of the normalized vector tanforming the orthonaoml basis determiend by PCA, oreder by decending variance of the data
%               The colomns of coeef are orthonormal {The columns are in the order of descending component variance, latent.}
%
%        scores is the representations of X in the principal component space. 
%               Each Row of score correspond to the vector of observation (the same veceor in the coresponding row of X) written in the new basis
%               Columns correspond to the components ("mixing" of observations
%
%               * We are intrested in the score (waight) of the most important dimension in each of vectors of Hot-Spots vlaues corresponding to a dimer (or a chain) 
%                 i.e. scores(:,1)
%              
%               PCA "centers" the data and uses the singular value decomposition (SVD) algorithm. 
%               => i.e. the mean scores of each feature (in our case hot spot) is substructed from all the scores of that feature  
%                  thus the 1st action of the PCA(mat_in) is: mat_in -> mat_in-mean(mat_in,1),  i.e. substract from each column its mean
%                  { Note that for the was Hot spot valuse are oreder in the main program and for SVD it shud be the mean of the rows i.e. mat_in-mean(mat_in,2)}
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, w1, w2]=residual_weighted_maps(map_in, sim_map_1, sim_map_2, mask_in)

if nargin==3
    mask_in=true(size(V));
end

V=map_in(mask_in);
u1=sim_map_1(mask_in);
u2=sim_map_2(mask_in);

w1 = - dot((V-u2),(u2-u1))/dot((u2-u1),(u2-u1));
w2= 1- w1;
res=V-w1*u1-w2*u2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [f1, loglike1, loglike2]=lk(sig, map_in, sim_map_1, sim_map_2, mask_in)

if nargin==3
    mask_in=true(size(V));
end

y=map_in(mask_in);
x1=sim_map_1(mask_in);
x2=sim_map_2(mask_in);
x3=(x1+x2)/2;


Const=dot(y-x1,y-x1)/(2*sig^2);

loglike1=0; % i.e. -dot(y-x1,y-x1)/(2*sig^2) + Const

loglike2=-dot(y-x2,y-x2)/(2*sig^2) + Const;
loglike3=-dot(y-x3,y-x3)/(2*sig^2) + Const;

%f1=exp(loglike3)/(exp(loglike1) + exp(loglike2) + exp(loglike3));

%f1=exp(loglike1)/(exp(loglike1) + exp(loglike2));
f1=1/(exp(loglike1) + exp(loglike2));


end

%% Scratch 

% 
% num_maps=size(dimers_list,1);
% for q1=1:num_maps
%     base_dir='C:\Users\frankg\OneDrive - BGU\Documents\_OutGoing\_Papers\BFR1-2\Comp\Chains\';
%     map_q1=ReadMRC([base_dir,dimers_list{q1,1}]);
%     if q1==1
%         sz_map=size(map_q1);
%         %dimers_TH    = nan(num_maps, 1);
%         dimers_notes = repmat({'empty'}      , num_maps, 1);
% 
%         dimers_TF    = repmat({false(sz_map)}, num_maps, 1);
%     end
%    
%     dimers_notes{q1} = dimers_list{q1,3};
%     dimers_TF{q1}    = map_q1 > dimers_list{q1,2};
% end

% num_maps=size(chains_list,1);
% for q1=1:num_maps
%     base_dir='C:\Users\frankg\OneDrive - BGU\Documents\_OutGoing\_Papers\BFR1-2\Comp\Chains\';
%     map_q1=ReadMRC([base_dir,chains_list{q1,1}]);
%     if q1==1
%         sz_map=size(map_q1);
%         %dimers_TH    = nan(num_maps, 1);
%         chains_notes = repmat({'empty'}      , num_maps, 1);
% 
%         chains_TF    = repmat({false(sz_map)}, num_maps, 1);
%     end
%    
%     chains_notes{q1} = chains_list{q1,3};
%     chains_TF{q1}    = map_q1 > chains_list{q1,2};
% end


% %% Building:
% % 1. PerDimer & PerChanin United HotSpot
% % 2. PerDimer & PerChain HS val Bfr1 hot-spots are postive and Bfr2 hot spots are negative 
% united_HS_dimers = repmat({false(sz_map)}, num_dimers, 1); % this will hold a TF map combaning all the hot spots in each dimer
% united_HS_chains = repmat({false(sz_map)}, num_dimers, 2); % this will hold a TF map combaning all the hot spots in each chain
% 
% HS_val_dimers_SymRlxMap=nan(num_dimers  , num_HS)';
% HS_val_chains_SymRlxMap=nan(num_dimers*2, num_HS)';
% dimers_order=repmat('z',num_dimers,1);
% 
% HS_val_dimers_ReSymMap=nan(num_dimers  , num_HS)';
% HS_val_chains_ReSymMap=nan(num_dimers*2, num_HS)';
% chains_order=repmat('z',num_dimers*2,1);
% 
% for q1=1:num_dimers % Iterating through all associated chain-pairs
%     united_HS_dimers_q1=false(sz_map);
%     united_HS_chain1_q1=false(sz_map);
%     united_HS_chain2_q1=false(sz_map);
% 
%     for q2=1:num_HS % iterating through all Hot Spots 
%         
%         % Finding the TF blob of hotspot q2 that belongs to the chain q1
%         [HS_q2_of_Bfr1_in_chain1_q1, vol_HS_q2_of_Bfr1_in_chain1_q1]=bigest_intersected_blob(chains1_TF{q1}, HS_Bfr1_TF{q2}); % Bfr1 hot-spot q2 in Chain1 in the dimer q1
%         [HS_q2_of_Bfr2_in_chain1_q1, vol_HS_q2_of_Bfr2_in_chain1_q1]=bigest_intersected_blob(chains1_TF{q1}, HS_Bfr2_TF{q2}); % Bfr2 hot-spot q2 in Chain1 in the dimer q1
% 
%         [HS_q2_of_Bfr1_in_chain2_q1, vol_HS_q2_of_Bfr1_in_chain2_q1]=bigest_intersected_blob(chains2_TF{q1}, HS_Bfr1_TF{q2}); % Bfr1 hot-spot q2 in Chain2 in the dimer q1
%         [HS_q2_of_Bfr2_in_chain2_q1, vol_HS_q2_of_Bfr2_in_chain2_q1]=bigest_intersected_blob(chains2_TF{q1}, HS_Bfr2_TF{q2}); % Bfr1 hot-spot q2 in Chain2 in the dimer q1
% 
%         % Consolidating the Hot-Spots TF maps of chains pairs composing a dimers 
%         HS_q2_of_Bfr1_in_dimer_q1 = HS_q2_of_Bfr1_in_chain1_q1 | HS_q2_of_Bfr1_in_chain2_q1;  % TF map of hot-spot q2 of Bfr1 in dimer q1 
%         vol_HS_q2_of_Bfr1_in_dimer_q1=sum(HS_q2_of_Bfr1_in_dimer_q1,"all");                   % "volume" i.e. the number of true voxels in the HS above 
% 
%         HS_q2_of_Bfr2_in_dimer_q1 = HS_q2_of_Bfr2_in_chain1_q1 | HS_q2_of_Bfr2_in_chain2_q1;  % TF map of hot-spot q2 of Bfr2 in dimer q1 
%         vol_HS_q2_of_Bfr2_in_dimer_q1=sum(HS_q2_of_Bfr2_in_dimer_q1,"all");
% 
%         comb_HS_q2_in_dimer_q1 = HS_q2_of_Bfr1_in_dimer_q1-HS_q2_of_Bfr2_in_dimer_q1;         % "Combaning" associated Bfr1 & Bfr2 the Hot-Spots of diemrs where Bfr1 is posstive and Bfr2 negative, overlaps are 0
%         vol_comb_HS_q2_in_dimer_q1 = sum(comb_HS_q2_in_dimer_q1~=0,"all");
%         united_HS_dimers_q1= united_HS_dimers_q1 | HS_q2_of_Bfr1_in_dimer_q1 | HS_q2_of_Bfr1_in_dimer_q1; % "Uniting" as above by | (or) operation 
%         
% 
%         comb_HS_q2_in_chain1_q1     = HS_q2_of_Bfr1_in_chain1_q1-HS_q2_of_Bfr2_in_chain1_q1; % "Combaning" associated Bfr1 & Bfr2 the Hot-Spots of chain1 where Bfr1 is posstive and Bfr2 negative, overlaps are 0
%         vol_comb_HS_q2_in_chain1_q1 = sum(comb_HS_q2_in_chain1_q1~=0,"all");
%         united_HS_chain1_q1= united_HS_chain1_q1 | HS_q2_of_Bfr1_in_chain1_q1 | HS_q2_of_Bfr2_in_chain1_q1; % "Uniting" as above by | (or) operation 
% 
% 
%         comb_HS_q2_in_chain2_q1     = HS_q2_of_Bfr1_in_chain2_q1-HS_q2_of_Bfr2_in_chain2_q1; % "Combaning" associated Bfr1 & Bfr2 the Hot-Spots of chain2 where Bfr1 is posstive and Bfr2 negative, overlaps are 0
%         vol_comb_HS_q2_in_chain2_q1 = sum(comb_HS_q2_in_chain2_q1~=0,"all");
%         united_HS_chain2_q1= united_HS_chain2_q1 | HS_q2_of_Bfr1_in_chain2_q1 | HS_q2_of_Bfr2_in_chain2_q1; % "Uniting" as above by | (or) operation 
% 
% 
%         % Each column is a dimer
%         % Each row    is a val of one of the hot-spots
%         HS_val_dimers_SymRlxMap(q2,q1)=sum(SymRlxMap.*comb_HS_q2_in_dimer_q1, "all")/vol_comb_HS_q2_in_dimer_q1;
%         HS_val_dimers_ReSymMap(q2,q1) =sum(ReSymMap .*comb_HS_q2_in_dimer_q1, "all")/vol_comb_HS_q2_in_dimer_q1;
% 
%         HS_val_dimers_SymRlxMap_NrmGlob(q2,q1)=sum(maps{1}.*comb_HS_q2_in_dimer_q1, "all")/vol_comb_HS_q2_in_dimer_q1;
%         HS_val_dimers_ReSymMap_NrmGlob(q2,q1) =sum(maps{2}.*comb_HS_q2_in_dimer_q1, "all")/vol_comb_HS_q2_in_dimer_q1;
% 
%         HS_val_dimers_SymRlxMap_NoNrm(q2,q1)=sum(MAPS{1}.*comb_HS_q2_in_dimer_q1, "all")/vol_comb_HS_q2_in_dimer_q1;
%         HS_val_dimers_ReSymMap_NoNrm(q2,q1) =sum(MAPS{2}.*comb_HS_q2_in_dimer_q1, "all")/vol_comb_HS_q2_in_dimer_q1;
% 
% 
%         
%         % Each column (q1 is constant) is a chain
%         % Each row    (q2 is constat)  is a val of one of the hot-spots
%         HS_val_chains_SymRlxMap(q2, 2*q1-1)=sum(SymRlxMap.*comb_HS_q2_in_chain1_q1, "all")/vol_comb_HS_q2_in_chain1_q1; 
%         HS_val_chains_SymRlxMap(q2, 2*q1  )=sum(SymRlxMap.*comb_HS_q2_in_chain2_q1, "all")/vol_comb_HS_q2_in_chain2_q1;
%         HS_val_chains_ReSymMap(q2, 2*q1-1 )=sum(ReSymMap .*comb_HS_q2_in_chain1_q1, "all")/vol_comb_HS_q2_in_chain1_q1;
%         HS_val_chains_ReSymMap(q2, 2*q1   )=sum(ReSymMap .*comb_HS_q2_in_chain2_q1, "all")/vol_comb_HS_q2_in_chain2_q1;
%         
%         chains_order(2*q1-1, 1) =dimers_couplings{q1,1};
%         chains_order(2*q1  , 1) =dimers_couplings{q1,2};
%     end % q2
%     united_HS_dimers{q1}   =united_HS_dimers_q1;
%     united_HS_chains{q1,1} =united_HS_chain1_q1;
%     united_HS_chains{q1,2} =united_HS_chain2_q1;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [U,S,V]=svd_proj_plot(HS_vals, HS_val_ReSym, dimers_names, ttl)
% %
% if nargin==3
%     ttl='';
% end
% 
% [U,S,V]=svd(HS_vals,"vector");
% ReSym=HS_val_ReSym;
% 
% proj_SV_1=(U(:,1)'*HS_vals)';
% h=plot_sorted_items(proj_SV_1, dimers_names, [ttl,' proj on 1st Singular Vector (explained variance)=', num2str(S(1).^2/sum(S.^2)*100)]);
% 
% 
% projReSym_SV_1=(U(:,1)'*ReSym)';
% %plot_sorted_items(projReSym_SV_1, dimers_names, [ttl,' Proj ReSym hot spots vals on 1st st Singular Vector ']);
% %mnMx=[min(projReSym_SV_1), max(projReSym_SV_1)]
% figure(h);
% hold on
% plot([0,13],[max(projReSym_SV_1), max(projReSym_SV_1)], '-.b')
% plot([0,13],[min(projReSym_SV_1), min(projReSym_SV_1)], '-.b')
% hold off
% end
