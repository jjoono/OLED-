%% Author : Jaehyeok Park (Eli)
% Transfer matrix simulator
% Date : 2016/07/11
% 2016 ProjectEli, All rights reserved.
% Inspired by Jinouk Song's TMF_whole function
% Date: 2016/08/15
% modified by JOSong from JHPark's work

%%
function TMF_birefringence_whole_p=TMF_birefringence_whole_p(no_bar,ne_bar,thickness,u,wavelength)

wavelength_num=length(wavelength);    % Number of wavelength point (default 301)

layer_num=size(no_bar,2);           % Number of layers

u_num=size(u,2);        % Number of u points

NL1=layer_num-1;

no_bar_u_num_1=repmat(no_bar,u_num,1);

cos_theta_p=sqrt(1-(repmat(reshape(u,wavelength_num*u_num,1),1,layer_num)./repmat(ne_bar,u_num,1)).^2);

% Physical branch: decay along +z requires imag(n*cos_theta)>=0 (see TMF_birefringence_whole.m).
flip_p=imag(repmat(ne_bar,u_num,1).*cos_theta_p)<0;
cos_theta_p(flip_p)=-cos_theta_p(flip_p);
clear flip_p

NL1vector=1:NL1;
NL1vector_plus1=NL1vector+1;

nj=no_bar_u_num_1(:,NL1vector);
nj1=no_bar_u_num_1(:,NL1vector_plus1);

cj_p=cos_theta_p(:,NL1vector);
cj1_p=cos_theta_p(:,NL1vector_plus1);

clear cos_theta_p

njcj_p=nj.*cj_p;
nj1cj=nj1.*cj_p;
njcj1=nj.*cj1_p;

clear nj nj1 cj_p cj1_p

denominator=nj1cj+njcj1;

l_r_p=(nj1cj-njcj1)./denominator;
l_t_p=2*njcj_p./denominator;

clear nj1cj njcj1 denominator

phasefactor_coefficient=2i*pi*repmat((1./wavelength)*thickness(NL1vector),u_num,1);

phasefactor_p=phasefactor_coefficient.*njcj_p;

clear njcj_p phasefactor_coefficient

minus_exp_p=exp(-phasefactor_p);
plus_exp_p=exp(phasefactor_p);

clear phasefactor_p

l_rt_p=l_r_p./l_t_p;

clear l_r_p

L_p=zeros(2,2,wavelength_num*u_num,NL1);

L_p(1,1,:,:)=minus_exp_p./l_t_p;
L_p(1,2,:,:)=minus_exp_p.*l_rt_p;
L_p(2,1,:,:)=plus_exp_p.*l_rt_p;
L_p(2,2,:,:)=plus_exp_p./l_t_p;

clear l_t_p l_rt_p minus_exp_p plus_exp_p

transfer_matrix=L_p(:,:,:,1);

for i=2:NL1
    
    transfer_matrix11=transfer_matrix(1,1,:).*L_p(1,1,:,i)+transfer_matrix(1,2,:).*L_p(2,1,:,i);
    transfer_matrix12=transfer_matrix(1,1,:).*L_p(1,2,:,i)+transfer_matrix(1,2,:).*L_p(2,2,:,i);
    transfer_matrix21=transfer_matrix(2,1,:).*L_p(1,1,:,i)+transfer_matrix(2,2,:).*L_p(2,1,:,i);
    transfer_matrix22=transfer_matrix(2,1,:).*L_p(1,2,:,i)+transfer_matrix(2,2,:).*L_p(2,2,:,i);
    
    transfer_matrix=horzcat(vertcat(transfer_matrix11,transfer_matrix21),vertcat(transfer_matrix12,transfer_matrix22));
    
end

clear L_p transfer_matrix11 transfer_matrix12 transfer_matrix21 transfer_matrix22

r_p=reshape(transfer_matrix(2,1,:)./transfer_matrix(1,1,:),wavelength_num,u_num);
t_p=reshape(1./transfer_matrix(1,1,:),wavelength_num,u_num);

clear transfer_matrix

TMF_birefringence_whole_p=struct('r_p',r_p,'t_p',t_p);