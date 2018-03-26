function d = kinectMesa(fvc, render_max)
  if nargin<2
    render_max = 0;
  end
  fvc1 = fvc;
  fvc.faces = reshape(1:size(fvc1.faces,1)*3, [size(fvc1.faces,2) size(fvc1.faces,1)])';
  t = reshape(fvc1.faces', [], 1);
  fvc.vertices = fvc1.vertices(t, :);

  %fvc.modelviewmatrix=create_lookat_matrix([-0.0244, -0.0259, 0], [-0.0244, -0.0259, 1], [0, 1, 0]);
  fvc.modelviewmatrix = [ 1.0000 0 0 0.0244; 0 1 0 0.0259; 0 0 1 0; 0 0 0 1];
  p = 7.8e-6;
  n = 519.46961112127485*p; f = 5;
  r = 1.0012*280*p; l = -1.0012*280*p; t = 213*p; b = -213*p;
  
  fvc.projectionmatrix=[2 * n / (r - l) 0 0 0; 0 2 * n / (t - b) 0 0; 0, 0, -(f + n) / (f - n), -(2 * f * n) / (f - n); 0 0 -1 0];
  fvc.vertices = double(fvc.vertices);
  fvc.viewport=[0 0 561 427];
  r=mesaPatch(fvc.projectionmatrix, fvc.modelviewmatrix, 561, 427, double(fvc.vertices'), uint32(fvc.faces'-1), uint32(render_max));
  r = flipud( r' ); d=zeros(427, 561); 
  if render_max
    m = (r~=0);
  else
    m = (r~=(2^32-1));
  end
  d(m) = 2*double(r(m))/(2^32-1)-1;
  d(m) = 2*n*f./(f+n-d(m)*(f-n));
end
