-- train script 

require 'sys'
require 'image'

local matio = require 'matio'

sampleSize = opt.batchSize
numberOfPasses = opt.numPasses

function getBatch_val(data, sampsize, count)
    -- select batch
    inputMat = torch.zeros(sampsize, data.inp:size(2), data.inp:size(3), data.inp:size(4))
    gtMat = torch.zeros(sampsize, data.gt:size(2), data.gt:size(3), data.gt:size(4))
    gt2Mat = torch.zeros(sampsize, data.gt2:size(2), data.gt2:size(3), data.gt2:size(4))
    
    for i = 1, sampsize do
        inputMat[{{i},{},{},{}}] = data.inp[{{count},{},{},{}}]
        gtMat[{{i},{},{}, {}}] = data.gt[{{count},{},{}, {}}]
        gt2Mat[{{i},{},{}, {}}] = data.gt2[{{count},{},{}, {}}]
        count = count + 1
    end

    return inputMat, gtMat, gt2Mat, count
end

function getBatch(data, sampsize, count, idx)
    -- select batch
    inputMat = torch.zeros(sampsize, data.inp:size(2), data.inp:size(3), data.inp:size(4))
    gtMat = torch.zeros(sampsize, data.gt:size(2), data.gt:size(3), data.gt:size(4))
    gtMsk = torch.zeros(sampsize, data.gt:size(2), data.gt:size(3), data.gt:size(4))
    
    gt2Mat = torch.zeros(sampsize, data.gt2:size(2), data.gt2:size(3), data.gt2:size(4))
    gt2Msk = torch.zeros(sampsize, data.gt2:size(2), data.gt2:size(3), data.gt2:size(4))

    for i = 1, sampsize do
        data_inp = data.inp[{{idx[count]},{},{},{}}]
        data_gt = data.gt[{{idx[count]},{},{}, {}}]
        data_gt2 = data.gt2[{{idx[count]},{},{}, {}}]
        -- data augmentation
        torch.seed() -- randomization
        -- flip
        local f_prob = torch.rand(1)
        if f_prob[1]>0.5 then
           data_inp = image.hflip(torch.reshape(data_inp, data.inp:size(2), data.inp:size(3), data.inp:size(4)))
           data_gt = image.hflip(torch.reshape(data_gt, data.gt:size(2), data.gt:size(3), data.gt:size(4)))  
          data_gt2 = image.hflip(torch.reshape(data_gt2, data.gt2:size(2), data.gt2:size(3), data.gt2:size(4))) 
          data_inp = torch.reshape(data_inp, 1, data.inp:size(2), data.inp:size(3), data.inp:size(4))
          data_gt = torch.reshape(data_gt, 1, data.gt:size(2), data.gt:size(3), data.gt:size(4))
	  data_gt2 = torch.reshape(data_gt2, 1, data.gt2:size(2), data.gt2:size(3), data.gt2:size(4))
        end
        -- gamma
        torch.seed() 
        local g_prob = torch.add(torch.mul(torch.rand(1),1.5), 0.5)
        data_inp = torch.pow(data_inp, g_prob[1])
        -- rotate
        torch.seed() 
        local r_prob = torch.add(torch.round(torch.mul(torch.rand(1),data.gt:size(4)-2)), 1) 
        data_inp = torch.cat(data_inp[{{},{},{},{r_prob[1]+1,data.inp:size(4)}}], data_inp[{{},{},{},{1,r_prob[1]}}], 4)
	data_gt = torch.cat(data_gt[{{},{},{},{r_prob[1]+1,data.gt:size(4)}}], data_gt[{{},{},{},{1,r_prob[1]}}], 4)
        data_gt2 = torch.cat(data_gt2[{{},{},{},{r_prob[1]+1,data.gt2:size(4)}}], data_gt2[{{},{},{},{1,r_prob[1]}}], 4)

        msk = data_gt:gt(0)
        msk2 = data_gt2:gt(0)
        inputMat[{{i},{},{},{}}] = data_inp
        gtMat[{{i},{},{}, {}}] = data_gt
        gtMsk[{{i},{},{}, {}}] = msk
        gt2Mat[{{i},{},{}, {}}] = data_gt2
        gt2Msk[{{i},{},{}, {}}] = msk2
        count = count + 1

        if count > tr_size then
            count = 1
            idx = torch.randperm(tr_size)
        end
    end

    return inputMat, gtMat, gtMsk, gt2Mat, gt2Msk, count, idx
end


function getValLoss()
    local valnumberOfPasses = torch.floor(pano_val.inp:size(1)/1)
    local loss = 0
    local valcount = 1
    --local out

    for i=1, valnumberOfPasses do

        --------------------- get mini-batch -----------------------
        inputMat, gtMat, gt2Mat, valcount = getBatch_val(pano_val, 1, valcount)
        ------------------------------------------------------------      

        -- forward
        inputMat = inputMat:cuda()
        gtMat = gtMat:cuda()
        gt2Mat = gt2Mat:cuda()
        --print('forward')
        output = model.core:forward(inputMat)
        --print(model.criterion:forward(output[1], gtMat))
        --print(model.criterion_2:forward(output[2], gt2Mat))

        loss = model.criterion:forward(output[1], gtMat) + model.criterion_2:forward(output[2], gt2Mat)+loss
        output = nil
        collectgarbage()
    end
    loss = loss / valnumberOfPasses

    return loss
end

-- do fwd/bwd and return loss, grad_params
function feval(x)
    if x ~= params then
        params:copy(x)
    end
    grad_params:zero()
    
    local loss = 0
    -- add for loop to increase mini-batch size
    for i=1, numberOfPasses do

        --------------------- get mini-batch -----------------------
        --inputMat, gtMat, gtMask = getBatch_rand(pano_tr, sampleSize)
        inputMat, gtMat, gtMsk, gt2Mat, gt2Msk, count, idx = getBatch(pano_tr, sampleSize, count, idx)
        ------------------------------------------------------------      

        -- forward
        inputMat = inputMat:cuda()
        gtMat = gtMat:cuda()
        gt2Mat = gt2Mat:cuda()

        output = model.core:forward(inputMat)
        
        --print(model.criterion:forward(output[1], gtMat))
        --print(model.criterion_2:forward(output[2], gt2Mat))
        
        loss = model.criterion:forward(output[1], gtMat) + model.criterion_2:forward(output[2], gt2Mat)+ loss
        
        -- backward
        loss_d_1 = model.criterion:backward(output[1], gtMat)
        loss_d_2 = model.criterion_2:backward(output[2], gt2Mat)
        
        gtMsk = torch.mul(gtMsk, 4)
        gtMsk = gtMsk:cuda()
        gtMsk_w = torch.cmul(loss_d_1, gtMsk)
        loss_d_1 = torch.add(gtMsk_w, loss_d_1)

        gt2Msk = torch.mul(gt2Msk, 4)
        gt2Msk = gt2Msk:cuda()
        gt2Msk_w = torch.cmul(loss_d_2, gt2Msk)
        loss_d_2 = torch.add(gt2Msk_w, loss_d_2)

        model.core:backward(inputMat, {loss_d_1, loss_d_2})
 
        output = nil
        loss_d_1 = nil
        loss_d_2 = nil
        gtMsk_w = nil
        gt2Msk_w = nil

        collectgarbage()
    end
    --print(loss_gt)
    
    grad_params:div(numberOfPasses)
    
    -- clip gradient element-wise
    grad_params:clamp(-10, 10)
    
    
    return loss/numberOfPasses, grad_params
end

losses = {} 
vallosses = {}
local optim_state = {opt.lr, opt.alpha, opt.epsilon}
local iterations = 8000
local minValLoss = 1/0
count = 1

idx = torch.randperm(pano_tr.inp:size(1))

for i = 1, iterations do
    
    
    model.core:training()
    local _, loss = optim.adam(feval, params, optim_state)
    --local _, loss = optim.rmsprop(feval, params, optim_state)

    print(string.format("update param, loss = %6.8f, gradnorm = %6.4e", loss[1], grad_params:clone():norm()))
    if i % 20 == 0 then
        print(string.format("iteration %4d, loss = %6.8f, gradnorm = %6.4e", i, loss[1], grad_params:norm()))
        model.core:evaluate()
        valLoss, output = getValLoss()
        vallosses[#vallosses + 1] = valLoss
        print(string.format("validation loss = %6.8f", valLoss))
        if minValLoss > valLoss then
            minValLoss = valLoss
            params_save = params:clone()
            nn.utils.recursiveType(params_save, 'torch.DoubleTensor')
            torch.save("./model/panofull_lay.t7", params_save:double())
            print("------- Model Saved --------")
        end
        losses[#losses + 1] = loss[1]
        
    end
end
