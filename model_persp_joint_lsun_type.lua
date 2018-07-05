
require 'torch'
require 'nn'
require 'nngraph'
require 'optim'
torch.setdefaulttensortype('torch.FloatTensor')
local model_utils=require 'model_utils'
require 'cunn'
require 'cudnn'

-- get training dataset
im_path = './data/LSUN_data/lsun_tr/img/'
ed_path = './data/LSUN_data/lsun_tr/edg/'
jc_path = './data/LSUN_data/lsun_tr/cor/'
jc_path2 = './data/LSUN_data/lsun_tr/cor_flip/'
id_path = './data/LSUN_data/lsun_tr/type/'
tr_name = torch.load('./data/LSUN_data/info_edg_stack_tr_lsun_640_d6_sig20_trname.t7')
pano_tr = {}
pano_tr.im_path = im_path
pano_tr.ed_path = ed_path
pano_tr.jc_path = jc_path
pano_tr.jc_path2 = jc_path2
pano_tr.tr_name = tr_name
pano_tr.id_path = id_path
tr_size = #tr_name
print(tr_size)

print('Uploaded training')

-- get validation dataset
img_val = torch.load('./data/LSUN_data/lsun_img_val.t7')
edg_val = torch.load('./data/LSUN_data/lsun_edg_val.t7')
junc_val = torch.load('./data/LSUN_data/lsun_juc_val.t7')
id_val = torch.load('./data/LSUN_data/lsun_id_val.t7')
pano_val = {}
pano_val.inp = img_val
pano_val.gt = edg_val
pano_val.gt2 = junc_val
pano_val.gt3 = id_val
val_size = pano_val.inp:size(1)
print(val_size)
print('Uploaded validation')

-- make model
model = {}

model.criterion = nn.BCECriterion():cuda()
model.criterion_2 = nn.BCECriterion():cuda()
model.criterion_3 = nn.CrossEntropyCriterion():cuda()

-- encoder
local input_x = nn.Identity()()
local conv1 = nn.SpatialConvolution(3,32,3,3,1,1,1,1)(input_x)
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
local deconv5_c = nn.SpatialConvolution(32*3,8,3,3,1,1,1,1)(unpool5_c)
local deconv6_sf_c = nn.Sigmoid()(deconv5_c)

-- refinement
local ref0 = nn.Reshape(2048*4*4)(pool7)
local ref1 = nn.Linear(2048*4*4, 1024)(ref0)
local ref1_relu = nn.ReLU(true)(ref1)
local ref2 = nn.Linear(1024, 256)(ref1_relu)
local ref2_relu = nn.ReLU(true)(ref2)
local ref3 = nn.Linear(256, 64)(ref2_relu)
local ref3_relu = nn.ReLU(true)(ref3)
local ref4 = nn.Linear(64, 11)(ref3_relu)

model.core = nn.gModule({input_x},{deconv6_sf, deconv6_sf_c, ref4})

model.core:cuda()

-- kaiming initialization
local method = 'kaiming'
model.core = require('weight-init')(model.core, method)

params, grad_params = model_utils.combine_all_parameters(model.core)

-- initialize from semantic segmentation
seg_params = torch.load('./model/perspfull_joint_lsun_pretrained.t7')
params[{{1, 94352587}}] = params[{{1, 94352587}}]:copy(seg_params)

print('start training')













