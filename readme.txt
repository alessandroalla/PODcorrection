 The code aims to approximate the following reaction-diffusion system

                   u_t = du \Delta u + f(u,v)
                   v_t = dv \Delta v + g(u,v)
          plus Boundary conditions and Initial conditions 
          (see eq (1) in https://arxiv.org/abs/2203.05998)
 using POD and POD-DEIM methods. We employ a correction term to stabilize
 the POD method and to improve the accuracy of the reduced method

 If you use this code please cite:
 Alessandro Alla, Angela Monti and Ivonne Sgura
 Adaptive POD-DEIM correction for Turing pattern approximation 
 in reaction-diffusion PDE systems
 accepted for Journal of Numerical Mathematics 
 (preprint https://arxiv.org/abs/2203.05998)


THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

The script can be run using the file main.m