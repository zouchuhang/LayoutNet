
require 'torch'
require 'nn'
require 'nngraph'
require 'optim'
torch.setdefaulttensortype('torch.FloatTensor')
local model_utils=require 'model_utils'
require 'cunn'
require 'cudnn'

-- get training dataset
print('loading training image ...')
img_tr = torch.load('./data/panoContext_img_train.t7')
print('loading training line map ...')
lne_tr = torch.load('./data/panoContext_line_train.t7')
print('loading training edge map ...')
edg_tr = torch.load('./data/panoContext_edge_train.t7')
print('loading training corner map ...')
junc_tr = torch.load('./data/panoContext_cor_train.t7')
print('loading training box param ...')
id_tr = torch.load('./data/panoContext_box_train.t7')
pano_tr = {}
pano_tr.inp = torch.cat(img_tr, lne_tr, 2)
pano_tr.gt = edg_tr
pano_tr.gt2 = junc_tr
pano_tr.gt3 = id_tr
tr_size = pano_tr.inp:size(1)
print(tr_size)
img_tr = nil
lne_tr = nil
edg_tr = nil
id_tr = nil
print('Uploaded training')

-- get validation dataset
img_val = torch.load('./data/panoContext_img_val.t7')
lne_val = torch.load('./data/panoContext_line_val.t7')
edg_val = torch.load('./data/panoContext_edge_val.t7')
junc_val = torch.load('./data/panoContext_cor_val.t7')
id_val = torch.load('./data/panoContext_box_val.t7')
pano_val = {}
pano_val.inp = torch.cat(img_val, lne_val, 2)
pano_val.gt = edg_val
pano_val.gt2 = junc_val
pano_val.gt3 = id_val
val_size = pano_val.inp:size(1)
print(val_size)
img_val = nil
lne_val = nil
edg_val = nil
id_val = nil
print('Uploaded validation')


-- make model
model = {}

model.criterion = nn.BCECriterion():cuda()
model.criterion_2 = nn.BCECriterion():cuda()
model.criterion_3 = nn.MSECriterion():cuda()

-- encoder
local input_x = nn.Identity()()
local conv1 = nn.SpatialConvolution(6,32,3,3,1,1,1,1)(input_x)
local conv1_relu = nn.ReLU(true)(conv1)
local pool1 = nn.SpatialMaxPooling(2,2,2,2)(conv1_relu)
local conv2 = nn.SpatialConvolution(32,64,3,3,1,1,1,1)(pool1)
local conv2_relu = nn.ReLU(true)(conv2)
local pool2 = nn.SpatialMaxPooling(2,2,2,2)(conv2_relu)
local conv3 = nn.SpatialConvolution(64,128,3,3,1,1,1,1)(pool2)
local conv3_relu = nn.ReLU(true)(conv3)
local pool3 = nn.SpatialMaxPooling(2,2,2,2)(conv3_relu)
local conv4 = nn.SpatialConvolution(128,256,3,3,1,1,1,1)(pool3)
local conv4_relu = nn.ReLU(true)(conv4)
local pool4 = nn.SpatialMaxPooling(2,2,2,2)(conv4_relu)
local conv5 = nn.SpatialConvolution(256,512,3,3,1,1,1,1)(pool4)
local conv5_relu = nn.ReLU(true)(conv5)
local pool5 = nn.SpatialMaxPooling(2,2,2,2)(conv5_relu)

local conv6 = nn.SpatialConvolution(512,1024,3,3,1,1,1,1)(pool5)
local conv6_relu = nn.ReLU(true)(conv6)
local pool6 = nn.SpatialMaxPooling(2,2,2,2)(conv6_relu)

local conv7 = nn.SpatialConvolution(1024,2048,3,3,1,1,1,1)(pool6)
local conv7_relu = nn.ReLU(true)(conv7)
local pool7 = nn.SpatialMaxPooling(2,2,2,2)(conv7_relu)

local unpool00 = nn.SpatialUpSamplingNearest(2)(pool7)
local deconv00 = nn.SpatialConvolution(2048,1024,3,3,1,1,1,1)(unpool00)
local deconv00_relu = nn.ReLU(true)(deconv00)

local unpool0_ = nn.JoinTable(2)({deconv00_relu, pool6})

local unpool0 = nn.SpatialUpSamplingNearest(2)(unpool0_)
local deconv0 = nn.SpatialConvolution(1024*2,512,3,3,1,1,1,1)(unpool0)
local deconv0_relu = nn.ReLU(true)(deconv0)

local unpool1_ = nn.JoinTable(2)({deconv0_relu, pool5})

local unpool1 = nn.SpatialUpSamplingNearest(2)(unpool1_)
local deconv1 = nn.SpatialConvolution(512*2,256,3,3,1,1,1,1)(unpool1)
local deconv1_relu = nn.ReLU(true)(deconv1)

local unpool2_ = nn.JoinTable(2)({deconv1_relu, pool4})

local unpool2 = nn.SpatialUpSamplingNearest(2)(unpool2_)
local deconv2 = nn.SpatialConvolution(256*2,128,3,3,1,1,1,1)(unpool2)
local deconv2_relu = nn.ReLU(true)(deconv2)

local unpool3_ = nn.JoinTable(2)({deconv2_relu, pool3})

local unpool3 = nn.SpatialUpSamplingNearest(2)(unpool3_)
local deconv3 = nn.SpatialConvolution(128*2,64,3,3,1,1,1,1)(unpool3)
local deconv3_relu = nn.ReLU(true)(deconv3)

local unpool4_ = nn.JoinTable(2)({deconv3_relu, pool2})

local unpool4 = nn.SpatialUpSamplingNearest(2)(unpool4_)
local deconv4 = nn.SpatialConvolution(64*2,32,3,3,1,1,1,1)(unpool4)
local deconv4_relu = nn.ReLU(true)(deconv4)

local unpool5_ = nn.JoinTable(2)({deconv4_relu, pool1})

local unpool5 = nn.SpatialUpSamplingNearest(2)(unpool5_)
local deconv5 = nn.SpatialConvolution(32*2,3,3,3,1,1,1,1)(unpool5)

local deconv6_sf = nn.Sigmoid()(deconv5)

-- joint part
local deconv00_c = nn.SpatialConvolution(2048,1024,3,3,1,1,1,1)(unpool00)
local deconv00_relu_c = nn.ReLU(true)(deconv00_c)
local unpool0_c = nn.JoinTable(2)({deconv00_relu_c, unpool0_})

local unpool0_c = nn.SpatialUpSamplingNearest(2)(unpool0_c)
local deconv0_c = nn.SpatialConvolution(1024*3,512,3,3,1,1,1,1)(unpool0_c)
local deconv0_relu_c = nn.ReLU(true)(deconv0_c)
local unpool1_c = nn.JoinTable(2)({deconv0_relu_c, unpool1_})

local unpool1_c = nn.SpatialUpSamplingNearest(2)(unpool1_c)
local deconv1_c = nn.SpatialConvolution(512*3,256,3,3,1,1,1,1)(unpool1_c)
local deconv1_relu_c = nn.ReLU(true)(deconv1_c)
local unpool2_c = nn.JoinTable(2)({deconv1_relu_c, unpool2_})

local unpool2_c = nn.SpatialUpSamplingNearest(2)(unpool2_c)
local deconv2_c = nn.SpatialConvolution(256*3,128,3,3,1,1,1,1)(unpool2_c)
local deconv2_relu_c = nn.ReLU(true)(deconv2_c)
local unpool3_c = nn.JoinTable(2)({deconv2_relu_c, unpool3_})

local unpool3_c = nn.SpatialUpSamplingNearest(2)(unpool3_c)
local deconv3_c = nn.SpatialConvolution(128*3,64,3,3,1,1,1,1)(unpool3_c)
local deconv3_relu_c = nn.ReLU(true)(deconv3_c)
local unpool4_c = nn.JoinTable(2)({deconv3_relu_c, unpool4_})

local unpool4_c = nn.SpatialUpSamplingNearest(2)(unpool4_c)
local deconv4_c = nn.SpatialConvolution(64*3,32,3,3,1,1,1,1)(unpool4_c)
local deconv4_relu_c = nn.ReLU(true)(deconv4_c)
local unpool5_c = nn.JoinTable(2)({deconv4_relu_c, unpool5_})

local unpool5_c = nn.SpatialUpSamplingNearest(2)(unpool5_c)
local deconv5_c = nn.SpatialConvolution(32*3,1,3,3,1,1,1,1)(unpool5_c)
local deconv6_sf_c = nn.Sigmoid()(deconv5_c)

-- 3D box
local deconv6_comb = nn.JoinTable(2)({deconv6_sf, deconv6_sf_c})
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

model.core = nn.gModule({input_x},{deconv6_sf, deconv6_sf_c, ref4})

model.core:cuda()

-- kaiming initialization
local method = 'kaiming'
model.core = require('weight-init')(model.core, method)

params, grad_params = model_utils.combine_all_parameters(model.core)

-- initialize from semantic segmentation
seg_params = torch.load('./model/panofull_lay_pretrained.t7')
params[{{1, 94347396}}] = params[{{1, 94347396}}]:copy(seg_params)

seg_params2 = torch.load('./model/panofull_box_pretrained.t7')
params[{{94347397, 94347396+18631262}}] = params[{{94347397, 94347396+18631262}}]:copy(seg_params2)
--print(seg_params[{{1,10}}])
--print(params[{{1,10}}])

print('start training')











