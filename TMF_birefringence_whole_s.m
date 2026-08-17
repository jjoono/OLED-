%% Author : Jaehyeok Park (Eli)
% Transfer matrix simulator
% Date : 2016/07/11
% 2016 ProjectEli, All rights reserved.
% Inspired by Jinouk Song's TMF_whole function
% Date: 2016/08/15
% modified by JOSong from JHPark's work

%%
function TMF_birefringence_whole_s=TMF_birefringence_whole_s(no_bar,ne_bar,thickness,u,wavelength)

wavelength_num=length(wavelength);    % Number of wavelength point (default 301)

layer_num=size(no_bar,2);           % Number of layers

u_num=size(u,2);        % Number of u points

NL1=layer_num-1;

no_bar_u_num_1=repmat(no_bar,u_num,1);

cos_theta_s=sqrt(1-(repmat(reshape(u,wavelength_num*u_num,1),1,layer_num)./no_bar_u_num_1).^2);

NL1vector=1:NL1;
NL1vector_plus1=NL1vector+1;

nj=no_bar_u_num_1(:,NL1vector);
nj1=no_bar_u_num_1(:,NL1vector_plus1);

cj_s=cos_theta_s(:,NL1vector);
cj1_s=cos_theta_s(:,NL1vector_plus1);

clear cos_theta_s

njcj_s=nj.*cj_s;
nj1cj1=nj1.*cj1_s;

clear nj nj1 cj_s cj1_s

l_r_s=(njcj_s-nj1cj1)./(njcj_s+nj1cj1);
l_t_s=1+l_r_s;

clear nj1cj1

phasefactor_coefficient=2i*pi*repmat((1./wavelength)*thickness(NL1vector),u_num,1);

phasefactor_s=phasefactor_coefficient.*njcj_s;

clear njcj_s phasefactor_coefficient

minus_exp_s=exp(-phasefactor_s);
plus_exp_s=exp(phasefactor_s);

clear phasefactor_s

l_rt_s=l_r_s./l_t_s;

clear l_r_s

L_s=zeros(2,2,wavelength_num*u_num,NL1);

L_s(1,1,:,:)=minus_exp_s./l_t_s;
L_s(1,2,:,:)=minus_exp_s.*l_rt_s;
L_s(2,1,:,:)=plus_exp_s.*l_rt_s;
L_s(2,2,:,:)=plus_exp_s./l_t_s;

clear l_t_s l_rt_s minus_exp_s plus_exp_s

transfer_matrix=L_s(:,:,:,1);

for i=2:NL1
    
    transfer_matrix11=transfer_matrix(1,1,:).*L_s(1,1,:,i)+transfer_matrix(1,2,:).*L_s(2,1,:,i);
    transfer_matrix12=transfer_matrix(1,1,:).*L_s(1,2,:,i)+transfer_matrix(1,2,:).*L_s(2,2,:,i);
    transfer_matrix21=transfer_matrix(2,1,:).*L_s(1,1,:,i)+transfer_matrix(2,2,:).*L_s(2,1,:,i);
    transfer_matrix22=transfer_matrix(2,1,:).*L_s(1,2,:,i)+transfer_matrix(2,2,:).*L_s(2,2,:,i);
    
    transfer_matrix=horzcat(vertcat(transfer_matrix11,transfer_matrix21),vertcat(transfer_matrix12,transfer_matrix22));
    
end

clear L_s transfer_matrix11 transfer_matrix12 transfer_matrix21 transfer_matrix22

r_s=reshape(transfer_matrix(2,1,:)./transfer_matrix(1,1,:),wavelength_num,u_num);
t_s=reshape(1./transfer_matrix(1,1,:),wavelength_num,u_num);

clear transfer_matrix

TMF_birefringence_whole_s=struct('r_s',r_s,'t_s',t_s);