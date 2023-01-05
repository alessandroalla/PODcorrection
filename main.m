% The code aims to approximate the following reaction-diffusion system
%
%                   u_t = du \Delta u + f(u,v)
%                   v_t = dv \Delta v + g(u,v)
%          plus Boundary conditions and Initial conditions 
%          (see eq (1) in https://arxiv.org/abs/2203.05998)
% using POD and POD-DEIM methods. We employ a correction term to stabilize
% the POD method and to improve the accuracy of the reduced method
%
% If you use this code please cite:
% Alessandro Alla, Angela Monti and Ivonne Sgura
% Adaptive POD-DEIM correction for Turing pattern approximation 
% in reaction-diffusion PDE systems
% accepted for Journal of Numerical Mathematics 
% (preprint https://arxiv.org/abs/2203.05998)
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
%FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
%COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
%IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
%CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

% Description of the variables
% model         choice of the model kinetics
%               1 --> FitzHugh-Nagumo
%               2 --> Schnakenberg
%               3 --> DIB
% choice        choice of the time interval
%               0 --> Time interval [0,T]
%               1 --> Zone I1 reactivity
%               2 --> Zone I2 stabilization
% data.U        POD basis for the unknown u
% data.V        POD basis for the unknown v
% R             dimension for the correction
% rdeim         number of interpolation points for DEIM
% r             dimension of the ROM
% ur, vr        reduced solutions found by ROMs

clc
clear all
close all

% Choose the model kinetics
model = 2;

% Compute snapshots
loading_snapshots

% Choose the time interval and set snapshots
choice = 1;
time_interval_choice

% Compute POD bases
PODbases

% Choose values for DEIM and correction
rankcorrection_choice

% Compute the correction term
[data.uR,data.vR,time_c] = pod_classic(model,R,data);

% Compute reduced solutions
for r = 1 : R
    % POD
    [ur,vr,time] = pod_classic(model,r,data); 
    u_pod = data.U(:,1:r)*ur(:,end); % projection back of reduced solution at final time T
    v_pod = data.V(:,1:r)*vr(:,end);
    ef_u_pod(r) = compute_error(u_pod,u_tf);
    ef_v_pod(r) = compute_error(v_pod,v_tf);
    time_pod(r) = time;
    
    % POD-DEIM
    [ur2,vr2,time2] = pod_deim_classic(model,r,rdeim,data); 
    u_deim = data.U(:,1:r)*ur2(:,end); % projection back of reduced solution at final time T
    v_deim = data.V(:,1:r)*vr2(:,end);
    ef_u_deim(r) = compute_error(u_deim,u_tf);
    ef_v_deim(r) = compute_error(v_deim,v_tf);
    time_deim(r) = time2; 
    
    % PODc
    [urc,vrc,time_corr] = pod_correction(model,r,R,data);
    u_podc = data.U(:,1:r)*urc(:,end); % projection back of reduced solution at final time T
    v_podc = data.V(:,1:r)*vrc(:,end);
    ef_u_podc(r) = compute_error(u_podc,u_tf);
    ef_v_podc(r) = compute_error(v_podc,v_tf);
    time_podc(r) = time_corr;
    
    %POD-DEIMc
    [urc2,vrc2,time2c] = pod_deim_correction(model,r,rdeim,R,data);
    u_deimc = data.U(:,1:r)*urc2(:,end); % projection back of reduced solution at final time T
    v_deimc = data.V(:,1:r)*vrc2(:,end);
    ef_u_deimc(r) = compute_error(u_deimc,u_tf);
    ef_v_deimc(r) = compute_error(v_deimc,v_tf);
    time_deimc(r) = time2c;
end    

% Generate plots
plots_generation
