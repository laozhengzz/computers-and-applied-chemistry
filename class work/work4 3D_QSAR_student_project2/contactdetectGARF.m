function [sum_com_dp_dp,sum_com_dp_indp,sum_com_indp_indp,sum_com_sc_sc] = contactdetectGARF(protein,ligand,Rcutoff,vdw_d)
rd = 0.005;
hnstate = 0.5/rd;
halva=hnstate;
halv=100;
halv1=10;
const=halva*2;
const1=halv1*2;
column_num=100;
Part_Matr_com=0;
Part_Matr_comA=0;

m=1;
n=1;
o=1;
p=1;

com_vdw_dp_dp=[];
com_vdw_dp_indp=[];
com_vdw_indp_indp=[];
com_vdw_sc_sc=[];


for i=1:1:size(protein,1)
    for j=1:1:size(ligand,1)
        if (sqrt((ligand(j,1)-protein(i,1))^2+(ligand(j,2)-protein(i,2))^2+(ligand(j,3)-protein(i,3))^2)<Rcutoff)
            
            a = protein(i,4);
            b = ligand(j,4);
            
            dist=norm(protein(i,1:3)-ligand(j,1:3));
            
            sd=[];
            sd=find(abs(vdw_d(:,1)-dist)<=rd/2);
            
            if isempty(sd)
                continue
            end
            
            order_num=(a-1)*34+b+1;
            
            com_vdw=vdw_d(sd,order_num);
            
            if ((protein(i,4) == 20 || (protein(i,4) >= 13 && protein(i,4) <= 16)) && (ligand(j,4) == 5 || (ligand(j,4) >= 7 && ligand(j,4) <= 13))) || ((ligand(j,4) == 5 || (ligand(j,4) >= 11 && ligand(j,4) <= 13)) && (protein(i,4) >= 13 && protein(i,4) <= 20))
                % dipole-dipole contacts
                com_vdw_dp_dp(m,1)=com_vdw;
                m=m+1;
            elseif ((protein(i,4) == 20 || (protein(i,4) >= 13 && protein(i,4) <= 16)) && ( (ligand(j,4) >= 2 && ligand(j,4) <= 4) || (ligand(j,4) >= 16 && ligand(j,4) <= 21) || (ligand(j,4) >= 29 && ligand(j,4) <= 34)) ) || ((ligand(j,4) == 5 || (ligand(j,4) >= 11 && ligand(j,4) <= 13)) && ( (protein(i,4) >= 1 && protein(i,4) <= 7) || (protein(i,4) >= 9 && protein(i,4) <= 10) || (protein(i,4) >= 21 && protein(i,4) <= 22)))
                % dipole-induceddipole contacts
                com_vdw_dp_indp(n,1)=com_vdw;
                n=n+1;
            elseif (( (protein(i,4) >= 1 && protein(i,4) <= 7) || (protein(i,4) >= 9 && protein(i,4) <= 10) || (protein(i,4) >= 21 && protein(i,4) <= 22)) && ( (ligand(j,4) >= 2 && ligand(j,4) <= 4) || (ligand(j,4) >= 16 && ligand(j,4) <= 21) || (ligand(j,4) >= 29 && ligand(j,4) <= 34)))
                % induceddipole-induceddipole contacts
                com_vdw_indp_indp(o,1)=com_vdw;
                o=o+1;
            elseif ( protein(i,4) == 8  && ( ligand(j,4) == 1 || (ligand(j,4) >= 26 && ligand(j,4) <= 28)) )
                % saturated carbon-saturated carbon contacts
                com_vdw_sc_sc(p,1)=com_vdw;
                p=p+1;
            end
        end
    end
end

sum_com_dp_dp=sum(com_vdw_dp_dp);

sum_com_dp_indp=sum(com_vdw_dp_indp);

sum_com_indp_indp=sum(com_vdw_indp_indp);

sum_com_sc_sc=sum(com_vdw_sc_sc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















%fclose(fp)



