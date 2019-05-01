function [sum_rep_dp_dp,sum_att_dp_dp,sum_rep_dp_indp,sum_att_dp_indp,sum_rep_indp_indp,sum_att_indp_indp,sum_rep_sc_sc,sum_att_sc_sc] = contactdetect(protein,ligand,Rcutoff)

m=1;
n=1;
o=1;
p=1;

contact_dp_dpA=[];
contact_dp_dpB=[];
dist_dp_dp=[];
rep_dp_dp=[];
att_dp_dp=[];
sum_rep_dp_dp=[];
sum_att_dp_dp=[];
contact_dp_indpA=[];
contact_dp_indpB=[];
dist_dp_indp=[];
rep_dp_indp=[];
att_dp_indp=[];
sum_rep_dp_indp=[];
sum_att_dp_indp=[];
contact_indp_indpA=[];
contact_indp_indpB=[];
dist_indp_indp=[];
rep_indp_indp=[];
att_indp_indp=[];
sum_rep_indp_indp=[];
sum_att_indp_indp=[];
contact_scA=[];
contact_scB=[];
dist_sc_sc=[];
rep_sc_sc=[];
att_sc_sc=[];
sum_rep_sc_sc=[];
sum_att_sc_sc=[];

for i=1:1:size(protein,1)
    for j=1:1:size(ligand,1)
        if (sqrt((ligand(j,1)-protein(i,1))^2+(ligand(j,2)-protein(i,2))^2+(ligand(j,3)-protein(i,3))^2)<Rcutoff)
            
            if ((protein(i,4) == 5 || (protein(i,4) >= 11 && protein(i,4) <= 13)) && (ligand(j,4) == 5 || (ligand(j,4) >= 7 && ligand(j,4) <= 13))) || ((ligand(j,4) == 5 || (ligand(j,4) >= 11 && ligand(j,4) <= 13)) && (protein(i,4) == 5 || (protein(i,4) >= 7 && protein(i,4) <= 13)))
                % dipole-dipole contacts
                contact_dp_dpA(m,:)=protein(i,:);
                contact_dp_dpB(m,:)=ligand(j,:);
                dist_dp_dp(m,1) = dist(contact_dp_dpA(m,:),contact_dp_dpB(m,:));
                rep_dp_dp(m,1) = (1/dist_dp_dp(m,1))^12;
                att_dp_dp(m,1) = (1/dist_dp_dp(m,1))^6;
                m=m+1;
            elseif ((protein(i,4) == 5 || (protein(i,4) >= 11 && protein(i,4) <= 13)) && ( (ligand(j,4) >= 2 && ligand(j,4) <= 4) || (ligand(j,4) >= 16 && ligand(j,4) <= 21) || (ligand(j,4) >= 29 && ligand(j,4) <= 34)) ) || ((ligand(j,4) == 5 || (ligand(j,4) >= 11 && ligand(j,4) <= 13)) && ( (protein(i,4) >= 2 && protein(i,4) <= 4) || (protein(i,4) >= 20 && protein(i,4) <= 21) || (protein(i,4) >= 29 && protein(i,4) <= 34)))
                % dipole-induceddipole contacts
                contact_dp_indpA(n,:)=protein(i,:);
                contact_dp_indpB(n,:)=ligand(j,:);
                dist_dp_indp(n,1) = dist(contact_dp_indpA(n,:),contact_dp_indpB(n,:));
                rep_dp_indp(n,1) = (1/dist_dp_indp(n,1))^12;
                att_dp_indp(n,1) = (1/dist_dp_indp(n,1))^6;
                n=n+1;
            elseif (( (protein(i,4) >= 2 && protein(i,4) <= 4) || (protein(i,4) >= 20 && protein(i,4) <= 21) || (protein(i,4) >= 29 && protein(i,4) <= 34)) && ( (ligand(j,4) >= 2 && ligand(j,4) <= 4) || (ligand(j,4) >= 16 && ligand(j,4) <= 21) || (ligand(j,4) >= 29 && ligand(j,4) <= 34)))
                % induceddipole-induceddipole contacts
                contact_indp_indpA(o,:)=protein(i,:);
                contact_indp_indpB(o,:)=ligand(j,:);
                dist_indp_indp(o,1) = dist(contact_indp_indpA(o,:),contact_indp_indpB(o,:));
                rep_indp_indp(o,1) = (1/dist_indp_indp(o,1))^12;
                att_indp_indp(o,1) = (1/dist_indp_indp(o,1))^6;
                o=o+1;
            elseif (( protein(i,4) == 1 || (protein(i,4) >= 26 && protein(i,4) <= 28)) && ( ligand(j,4) == 1 || (ligand(j,4) >= 26 && ligand(j,4) <= 28)) )
                % saturated carbon-saturated carbon contacts
                contact_scA(p,:)=protein(i,:);
                contact_scB(p,:)=ligand(j,:);
                dist_sc_sc(p,1) = dist(contact_scA(p,:),contact_scB(p,:));
                rep_sc_sc(p,1) = (1/dist_sc_sc(p,1))^12;
                att_sc_sc(p,1) = (1/dist_sc_sc(p,1))^6;
                p=p+1;
            end
        end
    end
end

sum_rep_dp_dp=sum(rep_dp_dp);
sum_att_dp_dp=sum(att_dp_dp);

sum_rep_dp_indp=sum(rep_dp_indp);
sum_att_dp_indp=sum(att_dp_indp);

sum_rep_indp_indp=sum(rep_indp_indp);
sum_att_indp_indp=sum(att_indp_indp);

sum_rep_sc_sc=sum(rep_sc_sc);
sum_att_sc_sc=sum(att_sc_sc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















%fclose(fp)



