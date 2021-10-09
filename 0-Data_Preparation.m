clc;
clear all;
clearvars;

stations={'Polmak'};
for jj=1:size(stations,2)
    clearvars -except jj stations
    filename = 'C:\Abolfazl\new2\PhD\Cases\Polmak.xlsx';
    [NUM,TXT,RAW] = xlsread(filename,stations{jj});
    tot_per= floor(size(NUM,1)/365); %total periods in one station
    intervals=[1,1.2,1.4,1.6,1.8,2];
    for j=1:size(intervals,2)
        k=1;
        for jjj=1:tot_per
            clear seq L_values seq_modified2 seq_modified seq_temp seq_temp2
            A(:,1) = NUM(k:jjj*365,1);
            A(:,2) = NUM(k:jjj*365,2);

            med = median(A(:,2));
            med=intervals(1,j)*med;
            miane(1:size(A,1),1)=med;

            %fig=plot(A(:,1),A(:,2),'-b');
            %hold 'on'

            %  plot(A(:,1),miane(:,1),'g');

            % consecutive events <= median value
            consenums = find(A(:,2)<=med);

            % plot(A(consenums,1),A(consenums,2),'.r')
            

            consenums=consenums' ;  %  array of consecutive numbers
            consenums(end+1)=2 ;  % adds new endpoint to very end of A so code picks up end of last group of consecutive values
            I_1=find(diff(consenums)~=1);  % finds where sequences of consecutive numbers end
            [m,n]=size(I_1);   % finds dimensions of I_1 i.e. how many sequences of consecutive numbers you have
            startpoint=1;    % sets start index at the first value in your array
            seq=cell(1,n) ; % had to preallocate because without, it only saved last iteration of the for loop below
                       % used n because this array is a row vector
                for i=1:n
                    End_Idx=I_1(i);   %set end index
                    seq{i}=consenums(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
                    startpoint=End_Idx+1;   %update start index for the next consecutive sequence
                end

            [row,col]=size(seq);
        
%                 for i2=1:col
%                     L_values(1,i2) = length(seq{i2});
%                 end
               %correction of the IIP by joining the close low discharge in
               %seq cell array
               
%                seq_temp=cell(1,col-1)
%                 if col>1
%                    for i4=1:col-1
%                        discharge_close=min(seq{i4+1})-max(seq{i4});
%                        dis_close(1,i4)=discharge_close;
%                    end
%                    seq_temp2=seq;
%                    [row2,col2]=size(dis_close);
%                    
%                    seq_modified=cell(1,col2);
                   %%%%%%%%%%%%%%%%%%%%%%%
%                    for i5=1:col2
%                        
%                        if dis_close(1,i5)<31 %25 days
%                           seq_modified{i5}= {[seq_temp2{i5},seq_temp2{i5+1}]};
%                        elseif dis_close(1,i5)>31 %25 days
%                            seq_modified{i5}= {[seq_temp2{i5+1}]};
%                        end
%                       
%                    end
                   %%%%%%%%%%%%%%%%%%%%%%%
%                    for i9=1:col2
%                        chk_var(1,i9)=isempty(seq_modified{1,i9});
%                    end
%                        if sum(chk_var)~=0
%                            for i8=1:col2
%                                L2=length(seq_temp2{i8});
%                            end
%                            [max_L2,Ind_max_L2]=max(L2);
%                            seq_modified={[seq_temp2{Ind_max_L2}]};
%                        end
% %                    [row3,col3]=size(seq_modified);
%                     [row5,col5]=size(seq_modified);
%                     if col5==1
%                         seq_modified0=seq_modified;
%                     else
%                         seq_modified0={[seq_modified{1,1:end}]};
%                     end
%                    [row4,col4]=size(seq_modified0{1,1});
%                    for i7=1:col4
%                    L1(i7)=length(seq_modified0{1,1}{i7});
%                    end
%                    [max_leng,ind_max_leng]=max(L1);
%                    
%                    seq_modified00={unique([seq_modified0{1,1}{ind_max_leng}])};
% %                    seq_modified000=unique(seq_modified00{1,1});
                       
                       
                
% %                     seq_modified2=seq_modified{1,1}(1,1);
%                     seq_modified2=seq_modified00;
%                     seq_modified2=seq_modified;
%                     [row,col]=size(seq_modified2);
%                     for i2=1:col
%                         L_values(1,i2) = length(seq_modified2{i2});
%                     end   
%                 else
                    seq_modified2=seq;
                    for i2=1:col
                        L_values(1,i2) = length(seq_modified2{i2});
                    end 
%                 end 
                   
               %two conditions for selecting the IIP one if the max discharge value
               %happens before 100 (example) and next is the highest
               %duration of IIP, if 
            duration = max(L_values); %days if ice influence period
                [max_dis_value,max_ind]=max(A(:,2));
                if max_ind<100
                    for i3=1:col
                        if length(seq_modified2{i3})==duration
                            ind=i3;
                            Ice_begin_1 = min(seq_modified2{ind});
                            Ice_end_1 = max(seq_modified2{ind});
                        end
                    end
                elseif max_ind >= 100
                    for i3=1:col
                        if length(seq_modified2{i3})==duration
                            ind=i3;
                            Ice_begin_1 = min(seq_modified2{ind});
                            Ice_end_1 = max(seq_modified2{ind});
                            if Ice_begin_1>max_ind
                                clear Ice_begin_1 Ice_end_1 row col i2 L_values
                                seq_modified2{i3}=[];
                                [row,col]=size(seq_modified2);
                                for i2=1:col
                                        L_values(1,i2) = length(seq_modified2{i2});
                                end
                                duration = max(L_values);
                                for i3=1:col
                                    if length(seq_modified2{i3})==duration
                                        ind=i3;
                                        Ice_begin_1 = min(seq_modified2{ind});
                                        Ice_end_1 = max(seq_modified2{ind});
                                    end
                                end
                            end
                        end
                    end
                    
                end
               
                
                %calculating the duration of ice influence period
                IIP_Duration(jjj,1)=Ice_end_1-Ice_begin_1+1;
                
                %getting the starting and ending point of the IIP(Ice
                %influence period)
                IIP_stdy(jjj,1)=Ice_begin_1;
                IIP_stdy(jjj,2)=Ice_end_1;
                if jjj>1
                    Ice_begin_1=Ice_begin_1+((jjj-1)*365);
                    Ice_end_1=Ice_end_1+((jjj-1)*365);
                end
                start_freezing=datestr(x2mdate(NUM(Ice_begin_1)));
                stop_freezing=datestr(x2mdate(NUM(Ice_end_1)));
                start_freezing=cellstr(start_freezing);
                stop_freezing=cellstr(stop_freezing);
                IIP_stdate(jjj,1)=start_freezing;
                IIP_stdate(jjj,2)=stop_freezing;
                
                %The first peak after IIP
                %if jjj=1%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if jjj>1
                    [M,I] = max(A((Ice_end_1-((jjj-1)*365))+1:365,2));
                    I_ind=I+(Ice_end_1-((jjj-1)*365));
                    Max_dis(jjj,2)=M;
                    Max_dis(jjj,1)=I_ind;
                    Max_dis(jjj,3)=I_ind-(Ice_end_1-((jjj-1)*365));
                else
                    [M,I] = max(A(Ice_end_1+1:365,2));
                    I_ind=I+Ice_end_1;
                    Max_dis(jjj,2)=M;
                    Max_dis(jjj,1)=I_ind;
                    Max_dis(jjj,3)=I_ind-Ice_end_1;
                end
                
                %
                k=k+365;
                inte={'10','12','14','16','18','20'};
                %temporary saving images to chck time series
                imagename=num2cell(1:85);
                cp={'r','k','c','m','y',[.8 .2 .6]};
                plot(A(:,1),A(:,2)); hold 'on'
                yline(med,'color',cp{j});
                if jjj==1
                    plot(A(Ice_begin_1:Ice_end_1,1),A(Ice_begin_1:Ice_end_1,2),'-r')
                end
                if jjj>1
                    plot(A(Ice_begin_1-((jjj-1)*365):Ice_end_1-((jjj-1)*365),1),A(Ice_begin_1-((jjj-1)*365):Ice_end_1-((jjj-1)*365),2),'-r')
                end
                filename2=num2str(jjj);
                saveas(gcf,filename2,'jpg')
                close
                
             
        end
        if jj==1 %this means the first station if jj==2 it is the second station
            
            xlswrite('C:\Abolfazl\Research\PhD\Cases\Kukko\startstopdate.xlsx',IIP_stdate,inte{j})
            xlswrite('C:\Abolfazl\Research\PhD\Cases\Kukko\day.xlsx',IIP_stdy,inte{j})
            xlswrite('C:\Abolfazl\Research\PhD\Cases\Kukko\maxdis.xlsx',Max_dis,inte{j})
            
        end
         if jj==1
%             
            xlswrite('C:\Abolfazl\Research\PhD\Cases\Kukko\IIP_Duration.xlsx',IIP_Duration,inte{j})
%            
        end

    end
end
