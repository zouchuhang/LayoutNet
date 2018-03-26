--

require 'torch'
require 'nn'
require 'nngraph'
require 'optim'
torch.setdefaulttensortype('torch.FloatTensor')
local model_utils=require 'model_utils'
require 'cunn'
require 'cudnn'

-- get training dataset
print('loading training corner map ...')
junc_tr = torch.load('./data/panoContext_cor_train.t7')
print('loading training edge map ...')
edg_tr = torch.load('./data/panoContext_edge_train.t7')
print('loading training box param ...')
id_tr = torch.load('./data/panoContext_box_train.t7')
pano_tr = {}
pano_tr.inp = junc_tr
pano_tr.inp2 = edg_tr
pano_tr.gt = id_tr
tr_size = pano_tr.inp:size(1)
print(tr_size)
junc_tr = nil
id_tr = nil
print('Uploaded training')

-- get validation dataset
junc_val = torch.load('./data/panoContext_cor_val.t7')
edg_val = torch.load('./data/panoContext_edge_val.t7')
id_val = torch.load('./data/panoContext_box_val.t7')
pano_val = {}
pano_val.inp = junc_val
pano_val.inp2 = edg_val
pano_val.gt = id_val
val_size = pano_val.inp:size(1)
print(val_size)

junc_val = nil
id_val = nil
print('Uploaded validation')


-- make model
model = {}

model.criterion = nn.MSECriterion():cuda()

-- refinement
local input_x = nn.Identity()()
local input_y = nn.Identity()()
local deconv6_comb = nn.JoinTable(2)({input_x, input_y})
local conv_r1 = nn.SpatialConvolution(4,8,3,3,1,1,1,1)(deconv6_comb)
local conv_r1_relu = nn.ReLU(true)(conv_r1)
local pool_r1 = nn.SpatialMaxPooling(2,2,2,2)(conv_r1_relu)
local conv_r2 = nn.SpatialConvolution(8,16,3,3,1,1,1,1)(pool_r1)
local conv_r2_relu = nn.ReLU(true)(conv_r2)
local pool_r2 = nn.SpatialMaxPooling(2,2,2,2)(conv_r2_relu)
local conv_r3 = nn.SpatialConvolution(16,32,3,3,1,1,1,1)(pool_r2)
local conv_r3_relu = nn.ReLU(true)(conv_r3)
local pool_r3 = nn.SpatialMaxPooling(2,2,2,2)(conv_r3_relu)
local conv_r4 = nn.SpatialConvolution(32,64,3,3,1,1,1,1)(pool_r3)
local conv_r4_relu = nn.ReLU(true)(conv_r4)
local pool_r4 = nn.SpatialMaxPooling(2,2,2,2)(conv_r4_relu)
local conv_r5 = nn.SpatialConvolution(64,128,3,3,1,1,1,1)(pool_r4)
local conv_r5_relu = nn.ReLU(true)(conv_r5)
local pool_r5 = nn.SpatialMaxPooling(2,2,2,2)(conv_r5_relu)
local conv_r6 = nn.SpatialConvolution(128,256,3,3,1,1,1,1)(pool_r5)
local conv_r6_relu = nn.ReLU(true)(conv_r6)
local pool_r6 = nn.SpatialMaxPooling(2,2,2,2)(conv_r6_relu)
local conv_r7 = nn.SpatialConvolution(256,512,3,3,1,1,1,1)(pool_r6)
local conv_r7_relu = nn.ReLU(true)(conv_r7)
local pool_r7 = nn.SpatialMaxPooling(2,2,2,2)(conv_r7_relu)

local ref0 = nn.Reshape(512*4*8)(pool_r7)
local ref1 = nn.Linear(512*4*8, 1024)(ref0)
local ref1_relu = nn.ReLU(true)(ref1)
local ref2 = nn.Linear(1024, 256)(ref1_relu)
local ref2_relu = nn.ReLU(true)(ref2)
local ref3 = nn.Linear(256, 64)(ref2_relu)
local ref3_relu = nn.ReLU(true)(ref3)
local ref4 = nn.Linear(64, 6)(ref3_relu)

model.core = nn.gModule({input_x, input_y},{ref4})

model.core:cuda()

-- kaiming initialization
local method = 'kaiming'
model.core = require('weight-init')(model.core, method)

params, grad_params = model_utils.combine_all_parameters(model.core)

print('start training')













