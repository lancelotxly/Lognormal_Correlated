clear all;clc;close all;
gamma_th=2;%0dBm
sigma_G=0.8;%lognormal��sigma����
PT=10:35;%���չ��ʣ�dB��
Pt=10.^(PT/10);%���չ��ʣ�W��
U_G=log(Pt.^0.5)-sigma_G^2;

P_EGC_Upper=[];
P_EGC_Lower=[];
sim_EGC=[]; 
arbtriary_upper_bound=[];
arbtriary_lower_bound=[];

for p=[0.1,0.8]
a=(1-sqrt(1-p^2))/p
U_X=U_G/(a+1);
sigma_X=sqrt(sigma_G^2/(a^2+1));

%Upper bound
OutageP_EGC_Upper=zeros(1,length(PT));%����һ���������ڴ���жϸ�������
parfor m=1:length(PT)%
    u_x=U_X(m);%����ѭ��ʹ�õ�mu_Gֵ
    OutageP_EGC_temp=qfunc(sqrt(2)*(u_x-1/(a+1)*log(sqrt(gamma_th/2)))/sigma_X);%���ۼ���ֵ
    OutageP_EGC_Upper(m)=OutageP_EGC_temp;%�����ۼ���ֵ��������
   
end
P_EGC_Upper=[P_EGC_Upper;OutageP_EGC_Upper];

%the lower bound
OutageP_EGC_Lower=zeros(1,length(PT));%����һ���������ڴ���жϸ�������
parfor m=1:length(PT)%
    u_x=U_X(m);%����ѭ��ʹ�õ�mu_Gֵ
    Q=qfunc((u_x-1/(a+1)*log(sqrt(gamma_th/2)))/sigma_X);
    OutageP_EGC_temp=Q^2;%���ۼ���ֵ
    OutageP_EGC_Lower(m)=OutageP_EGC_temp;%�����ۼ���ֵ��������
end
P_EGC_Lower=[P_EGC_Lower;OutageP_EGC_Lower];

OutageP_sim_EGC=0;%Simulation
M=100;N=500;%�������������
for m=1:M
m%��ʾ��ǰѭ������
OutageP_sim_vec_EGC_temp=zeros(1,length(PT));%����һ���������ڴ�ŷ���õ����жϸ�������
parfor k=1:length(PT)
    u_x=U_X(k);    %����ѭ��ʹ�õ�mu_Gֵ
    %G_vec=(sigma_G*randn(2,N)+mu_G);  %����L*N��Gauss�������
    X_1=sigma_X*randn(1,N)+u_x;
    X_2=sigma_X*randn(1,N)+u_x;
    G_vec_1=a*X_1+X_2;
    G_vec_2=X_1+a*X_2;
    G_vec=[G_vec_1;G_vec_2];
    OutageP_sim_temp=sum(double(sum(exp(G_vec))<(2*gamma_th)^0.5))/N;%������Щ����������жϸ���
    OutageP_sim_vec_EGC_temp(k)=OutageP_sim_temp;%������õ��Ľ����������
end
OutageP_sim_EGC=OutageP_sim_EGC+OutageP_sim_vec_EGC_temp;%
end
OutageP_sim_EGC=OutageP_sim_EGC/M;%��һ���ǶԶ�η���õ����жϸ�������ȡƽ������Ϊ�ڴ治���Դ�M*N�������

sim_EGC=[sim_EGC;OutageP_sim_EGC];

%arbtriary upper-bound
OutageP_EGC_Upper_X=zeros(1,length(PT));%����һ���������ڴ���жϸ�������
parfor m=1:length(PT)%
    u_x=U_X(m);%����ѭ��ʹ�õ�mu_Gֵ
    if p==0.1
        x_m1=-1.508;
        x_m2=0.6535;
    else if p==0.8
            x_m1=1.3;
            x_m2=-1.671;
        end
    end
    tempA_X=abs(u_x-x_m2+(a*exp(a*x_m1+x_m2)+exp(x_m1+a*x_m2))/(exp(a*x_m1+x_m2)+a*exp(x_m1+a*x_m2))*(u_x-x_m1));
    tempB_X=sqrt(sigma_X^2*(1+((a*exp(a*x_m1+x_m2)+exp(x_m1+a*x_m2))/(exp(a*x_m1+x_m2)+a*exp(x_m1+a*x_m2)))^2));
    OutageP_EGC_temp=qfunc(tempA_X/tempB_X)%���ۼ���ֵ
    OutageP_EGC_Upper_X(m)=OutageP_EGC_temp;%�����ۼ���ֵ��������
end
arbtriary_upper_bound=[arbtriary_upper_bound;OutageP_EGC_Upper_X];

%arbtriary lower-bound
OutageP_EGC_Lower_X=zeros(1,length(PT));%����һ���������ڴ���жϸ�������
parfor m=1:length(PT)%
    u_x=U_X(m);%����ѭ��ʹ�õ�mu_Gֵ
  if p==0.1
        x_m1=-1.508;
        x_m2=0.6535;
    else if p==0.8
            x_m1=1.3;
            x_m2=-1.671;
        end
    end
    Q_m1=qfunc((u_x-x_m1)/sigma_X);
    Q_m2=qfunc((u_x-x_m2)/sigma_X);
    OutageP_EGC_temp=Q_m1*Q_m2;%���ۼ���ֵ
    OutageP_EGC_Lower_X(m)=OutageP_EGC_temp;%�����ۼ���ֵ��������
end
arbtriary_lower_bound=[arbtriary_lower_bound;OutageP_EGC_Lower_X];
end

%%��ͼ
H1_1=semilogy(PT,P_EGC_Upper(1,:),'b--','linewidth',1); hold on;%��ͼ
H1_2=semilogy(PT,P_EGC_Upper(2,:),'r--','linewidth',1); hold on;%��ͼ
%H1_3=semilogy(PT,P_EGC_Upper(3,:),'k--','linewidth',1); hold on;%��ͼ
% 
H2_1=semilogy(PT,P_EGC_Lower(1,:),'b--','linewidth',1); hold on;%��ͼ
H2_2=semilogy(PT,P_EGC_Lower(2,:),'r--','linewidth',1); hold on;%��ͼ
%H2_3=semilogy(PT,P_EGC_Lower(3,:),'k--','linewidth',1); hold on;%��ͼ

H3_1=semilogy(PT,sim_EGC(1,:),'bd','linewidth',1); hold on;
H3_2=semilogy(PT,sim_EGC(2,:),'rs','linewidth',1); hold on;
%H3_3=semilogy(PT,sim_EGC(3,:),'ko','linewidth',1); hold on;
% 
H4_1=semilogy(PT,arbtriary_upper_bound(1,:),'b-','linewidth',1); hold on;%��ͼ
H4_2=semilogy(PT,arbtriary_upper_bound(2,:),'r-','linewidth',1); hold on;%��ͼ
% H4_3=semilogy(PT,arbtriary_upper_bound(3,:),'k-','linewidth',1); hold on;%��ͼ
% %  
H5_1=semilogy(PT,arbtriary_lower_bound(1,:),'b-','linewidth',1); hold on;%��ͼ
H5_2=semilogy(PT,arbtriary_lower_bound(2,:),'r-','linewidth',1); hold on;%��ͼ
% %H5_3=semilogy(PT,arbtriary_lower_bound(3,:),'k-','linewidth',1); hold on;%��ͼ
% % %  
grid on;
grid minor;

legend([H1_1,H1_2,H4_1,H4_2,H3_1,H3_2],'P_{{upper},{{x}_{m}^{*}}}^{EGC} or P_{{lower},{x_{m}^*}}^{EGC},\rho=0.1',...
'P_{{upper},{x}_{m}^{*}}^{EGC} or P_{{lower},{x}_{m}^{*}}^{EGC}, \rho=0.8','P_{{upper},{x_m}}^{EGC} or P_{{lower},{x_m}}^{EGC}, \rho=0.1',...
'P_{{upper},{x_m}}^{EGC} or P_{{lower},{x_m}}^{EGC}, \rho=0.8','Sim. P_{out}^{EGC}, \rho=0.1','Sim. P_{out}^{EGC}, \rho=0.8',...
 'Location','southwest');
  
% legend([H1_1,H1_2,H1_3,H3_1,H3_2,H3_3],'P_{Upper}^{EGC} or P_{Lower}^{EGC}, \rho=0.01','P_{Upper}^{EGC} or P_{Lower}^{EGC}, \rho=0.6',...
% 'P_{Upper}^{EGC} or P_{Lower}^{EGC}, \rho=0.99','Sim. P_{out}^{EGC}, \rho=0.01','Sim. P_{out}^{EGC}, \rho=0.6','Sim. P_{out}^{EGC},  \rho=0.99',...
%  'Location','southwest');

xlabel({'$\overline{E} (dB)$'},'Interpreter','latex');
ylabel('Outage Probability');




