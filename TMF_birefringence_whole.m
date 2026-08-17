%% Author : Jaehyeok Park (Eli)
% Transfer matrix simulator
% Date : 2016/07/11
% 2016 ProjectEli, All rights reserved.
% Inspired by Jinouk Song's TMF_whole function
% Date: 2016/08/15
% modified by JOSong from JHPark's work

%%
function TMF_birefringence_whole=TMF_birefringence_whole(no_bar,ne_bar,thickness,u,wavelength)

wavelength_num=length(wavelength);    % Number of wavelength point (default 301)

layer_num=size(no_bar,2);           % Number of layers

u_num=size(u,2);        % Number of u points

NL1=layer_num-1;

no_bar_u_num_1=repmat(no_bar,u_num,1);
ne_bar_u_num_1=repmat(ne_bar,u_num,1);

cos_theta_p=sqrt(1-(repmat(reshape(ne_bar(:,1)*u,wavelength_num*u_num,1),1,layer_num)./ne_bar_u_num_1).^2);
cos_theta_s=sqrt(1-(repmat(reshape(no_bar(:,1)*u,wavelength_num*u_num,1),1,layer_num)./no_bar_u_num_1).^2);

NL1vector=1:NL1;
NL1vector_plus1=NL1vector+1;

nj=no_bar_u_num_1(:,NL1vector);
nj1=no_bar_u_num_1(:,NL1vector_plus1);

cj_p=cos_theta_p(:,NL1vector);
cj1_p=cos_theta_p(:,NL1vector_plus1);

% clear cos_theta_p

cj_s=cos_theta_s(:,NL1vector);
cj1_s=cos_theta_s(:,NL1vector_plus1);

% clear cos_theta_s

njcj_p=nj.*cj_p;
njcj_s=nj.*cj_s;
nj1cj=nj1.*cj_p;
njcj1=nj.*cj1_p;
nj1cj1=nj1.*cj1_s;

clear nj nj1 cj_p cj1_p cj_s cj1_s

denominator=nj1cj+njcj1;

l_r_p=(nj1cj-njcj1)./denominator;
l_t_p=2*njcj_p./denominator;

clear nj1cj njcj1 denominator

l_r_s=(njcj_s-nj1cj1)./(njcj_s+nj1cj1);
l_t_s=1+l_r_s;

clear nj1cj1

phasefactor_coefficient=2i*pi*repmat((1./wavelength)*thickness(NL1vector),u_num,1);

phasefactor_p=phasefactor_coefficient.*njcj_p;
phasefactor_s=phasefactor_coefficient.*njcj_s;

clear njcj_p njcj_s phasefactor_coefficient

minus_exp_p=exp(-phasefactor_p);
plus_exp_p=exp(phasefactor_p);

clear phasefactor_p

minus_exp_s=exp(-phasefactor_s);
plus_exp_s=exp(phasefactor_s);

clear phasefactor_s

l_rt_p=l_r_p./l_t_p;
l_rt_s=l_r_s./l_t_s;

clear l_r_p l_r_s

L_p=zeros(2,2,wavelength_num*u_num,NL1);
L_s=zeros(2,2,wavelength_num*u_num,NL1);

L_p(1,1,:,:)=minus_exp_p./l_t_p;
L_p(1,2,:,:)=minus_exp_p.*l_rt_p;
L_p(2,1,:,:)=plus_exp_p.*l_rt_p;
L_p(2,2,:,:)=plus_exp_p./l_t_p;

clear l_t_p l_rt_p minus_exp_p plus_exp_p

L_s(1,1,:,:)=minus_exp_s./l_t_s;
L_s(1,2,:,:)=minus_exp_s.*l_rt_s;
L_s(2,1,:,:)=plus_exp_s.*l_rt_s;
L_s(2,2,:,:)=plus_exp_s./l_t_s;

clear l_t_s l_rt_s minus_exp_s plus_exp_s

L_ps=cat(3,L_p,L_s);

clear L_p L_s

transfer_matrix=L_ps(:,:,:,1);

for i=2:NL1
    
    transfer_matrix11=transfer_matrix(1,1,:).*L_ps(1,1,:,i)+transfer_matrix(1,2,:).*L_ps(2,1,:,i);
    transfer_matrix12=transfer_matrix(1,1,:).*L_ps(1,2,:,i)+transfer_matrix(1,2,:).*L_ps(2,2,:,i);
    transfer_matrix21=transfer_matrix(2,1,:).*L_ps(1,1,:,i)+transfer_matrix(2,2,:).*L_ps(2,1,:,i);
    transfer_matrix22=transfer_matrix(2,1,:).*L_ps(1,2,:,i)+transfer_matrix(2,2,:).*L_ps(2,2,:,i);
    
    transfer_matrix=horzcat(vertcat(transfer_matrix11,transfer_matrix21),vertcat(transfer_matrix12,transfer_matrix22));
    
end

% clear L_ps transfer_matrix11 transfer_matrix12 transfer_matrix21 transfer_matrix22

r=transfer_matrix(2,1,:)./transfer_matrix(1,1,:);
t=1./transfer_matrix(1,1,:);

% clear transfer_matrix

r_p=reshape(r(1:wavelength_num*u_num),wavelength_num,u_num);
r_s=reshape(r(1+wavelength_num*u_num:2*wavelength_num*u_num),wavelength_num,u_num);

clear r

t_p=reshape(t(1:wavelength_num*u_num),wavelength_num,u_num);
t_s=reshape(t(1+wavelength_num*u_num:2*wavelength_num*u_num),wavelength_num,u_num);

clear t

% TMF_birefringence_whole=struct('r_p',r_p,'t_p',t_p,'r_s',r_s,'t_s',t_s, 'transfer_matrix11', transfer_matrix, 'transfer_matrix12' ,transfer_matrix12, 'transfer_matrix21', transfer_matrix21, 'transfer_matrix22', transfer_matrix22);
TMF_birefringence_whole=struct('r_p',r_p,'t_p',t_p,'r_s',r_s,'t_s',t_s, 'transfer_matrix', transfer_matrix, 'L_ps', L_ps, 'cos_p', cos_theta_p, 'cos_s', cos_theta_s);