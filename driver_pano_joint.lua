-- driver

local cmd = torch.CmdLine()
cmd:text()
cmd:text('Script for training LayoutNet corner & edge branch.')

cmd:option('-lr' , 1e-4, 'learning rate')
cmd:option('-alpha' , 0.95, 'adam, alpha')
cmd:option('-epsilon' , 1e-6, 'adam, epsilon')
cmd:option('-batchSize' , 1, 'mini batch size')
cmd:option('-numPasses' , 5, 'number of passes, works as batch size to avoid out of memory')
--cmd:option('-valData' , './full_pano_processed/pano_val.t7', 'filepath for validation data')
--cmd:option('-trainData' , './full_pano_processed/pano_tr.t7', 'filepath for training data')

cmd:text()
opt = cmd:parse(arg)

dofile('model_pano_joint.lua')
dofile('train_pano_joint.lua')
