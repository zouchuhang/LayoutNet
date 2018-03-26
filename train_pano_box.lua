-- train script

require 'sys'
require 'image'

local matio = require 'matio'

sampleSize = opt.batchSize
numberOfPasses = opt.numPasses

function getBatch_val(data, sampsize, count)
    -- select batch
    inputMat = torch.zeros(sampsize, data.inp:size(2), data.inp:size(3), data.inp:size(4))
    inputMat2 = torch.zeros(sampsize, data.inp2:size(2), data.inp2:size(3), data.inp2:size(4))
    gtMat = torch.zeros(sampsize, 6)

    for i = 1, sampsize do
        inputMat[{{i},{},{},{}}] = data.inp[{{count},{},{},{}}]
        inputMat2[{{i},{},{},{}}] = data.inp2[{{count},{},{},{}}]
        gtMat[{{i},{}}] = data.gt[{{count},{1}, {1,6}}]
        count = count + 1
    end

    return inputMat, inputMat2, gtMat, count
end

function getBatch(data, sampsize, count, idx)
    -- select batch
    inputMat = torch.zeros(sampsize, data.inp:size(2), data.inp:size(3), data.inp:size(4))
    inputMat2 = torch.zeros(sampsize, data.inp2:size(2), data.inp2:size(3), data.inp2:size(4))
    gtMat = torch.zeros(sampsize,6)

    for i = 1, sampsize do
        data_inp = data.inp[{{idx[count]},{},{},{}}]
        data_inp2 = data.inp2[{{idx[count]},{},{},{}}]
        data_gt = data.gt[{{idx[count]},{1}, {1,6}}]
        -- data augmentation
        torch.seed() -- randomization
        -- rotate
        torch.seed()
        local r_prob = torch.add(torch.round(torch.mul(torch.rand(1),data.gt:size(2)-2)), 1)
        if r_prob[1] > 0 then 
        data_inp = torch.cat(data_inp[{{},{},{},{data.inp:size(4) - r_prob[1]+1,data.inp:size(4)}}], data_inp[{{},{},{},{1,data.inp:size(4) - r_prob[1]}}], 4)
        data_inp2 = torch.cat(data_inp2[{{},{},{},{data.inp2:size(4) - r_prob[1]+1,data.inp2:size(4)}}], data_inp2[{{},{},{},{1,data.inp2:size(4) - r_prob[1]}}], 4)
        data_gt = data.gt[{{idx[count]},{r_prob[1]+1}, {1,6}}]
        end
        -- flip
        local f_prob = torch.rand(1)
        if f_prob[1]>0.5 then
           data_inp = image.hflip(torch.reshape(data_inp, data.inp:size(2), data.inp:size(3), data.inp:size(4))) 
          data_inp2 = image.hflip(torch.reshape(data_inp2, data.inp2:size(2), data.inp2:size(3), data.inp2:size(4)))
          data_inp = torch.reshape(data_inp, 1, data.inp:size(2), data.inp:size(3), data.inp:size(4))
          data_inp2 = torch.reshape(data_inp2, 1, data.inp2:size(2), data.inp2:size(3), data.inp2:size(4))
          data_gt = data.gt[{{idx[count]},{r_prob[1]+1}, {7,12}}]
        end
        -- gamma
        --torch.seed() 
        --local g_prob = torch.add(torch.mul(torch.rand(1),1.5), 0.5)
        --data_inp = torch.pow(data_inp, g_prob[1])

        inputMat[{{i},{},{},{}}] = data_inp
        inputMat2[{{i},{},{},{}}] = data_inp2
 
        gtMat[{{i},{}}] = data_gt
        count = count + 1

        if count > tr_size then
            count = 1
            idx = torch.randperm(tr_size)
        end
    end

    return inputMat, inputMat2, gtMat, count, idx
end


function getValLoss()
    local valnumberOfPasses = torch.floor(pano_val.inp:size(1)/1)
    local loss = 0
    local valcount = 1
    --local out

    for i=1, valnumberOfPasses do

        --------------------- get mini-batch -----------------------
        inputMat, inputMat2, gtMat, valcount = getBatch_val(pano_val, 1, valcount)
        ------------------------------------------------------------      

        -- forward
        inputMat = inputMat:cuda()
        inputMat2 = inputMat2:cuda()
        gtMat = gtMat:cuda()
      
        --print('forward')
        output = model.core:forward({inputMat, inputMat2})

        loss = model.criterion:forward(output, gtMat) + loss
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
        inputMat, inputMat2, gtMat, count, idx = getBatch(pano_tr, sampleSize, count, idx)
        ------------------------------------------------------------      

        -- forward
        inputMat = inputMat:cuda()
        inputMat2 = inputMat2:cuda()
        gtMat = gtMat:cuda()

        output = model.core:forward({inputMat, inputMat2})
        
        loss = model.criterion:forward(output, gtMat) + loss

        -- backward
        loss_d_1 = model.criterion:backward(output, gtMat)

        model.core:backward({inputMat,inputMat2}, loss_d_1)
 
        output = nil
        loss_d_1 = nil

        collectgarbage()
    end
    
    grad_params:div(numberOfPasses)
    
    -- clip gradient element-wise
    grad_params:clamp(-10, 10)
    
    
    return loss/numberOfPasses, grad_params
end

losses = {} 
vallosses = {}
local optim_state = {opt.lr, opt.alpha, opt.epsilon}

local iterations = 9000
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
            torch.save("./model/panofull_box.t7", params_save:double())
            print("------- Model Saved --------")
        end
        losses[#losses + 1] = loss[1]
    end
end
