nx = 32;
ny = 32;

N = 5000;
load ('mnist.mat');
XTrain = training.images(:,:,1:N);
YTrain = training.labels(1:N) + 1;

files = dir('data/5000_mnist/*.mat');
total_files = length(files);
X = zeros(nx*nx + 28*28, total_files);
Y = zeros(1, total_files);

stream = RandStream('mt19937ar','seed', 1);
RandStream.setGlobalStream(stream);
shuffle_idxs = randperm(total_files);
for i=1:total_files
    f = files(i);
    idx = sscanf(f.name, "phases_%d.mat");
    load(['data/5000_mnist/' f.name])
    inputs = 1 - 2 * XTrain(:,:,idx);
    X(:,shuffle_idxs(i)) = [reshape(phases, nx*ny, 1); reshape(inputs, 28*28, 1)];
    Y(shuffle_idxs(i)) = YTrain(idx);
end
test_classes = 1:10; %[4,9];
nclasses = length(test_classes);

YY = Y; %(Y==test_classes(1) | Y==test_classes(2));
XX = X; %(:, Y==test_classes(1) | Y==test_classes(2));

%ntrain = round(nfiles * 0.8);   
%ntest = nfiles - ntrain;
ntrain = round(length(YY)*0.8);
ntest = length(YY) - ntrain;

train_idxs = 1:ntrain;
test_idxs = ntrain+1:length(YY); 

alpha = 1e3;
noutputs = size(XX,1);
wfrac = 1;
nweights = round(noutputs*wfrac);
widx = 1:nweights; %randperm(noutputs, nweights);

% Train output weights
Ytrain = Y(train_idxs);
votes = zeros(ntrain, nclasses);
W = zeros(nweights, nclasses, nclasses, 2);
for c1=1:nclasses-1
    for c2=c1+1:nclasses
        img_idxs = Ytrain==c1 | Ytrain==c2;
        Yc1c2 = Ytrain(img_idxs);
        Xc1c2 = XX(widx, img_idxs);
        Yc1 = (Yc1c2==c1)*1 - (Yc1c2==c2)*1;
        Yc2 = (Yc1c2==c2)*1 - (Yc1c2==c1)*1;
        Wc1 = Yc1 * Xc1c2' / (Xc1c2*Xc1c2' + alpha*eye(nweights,nweights));
        Wc2 = Yc2 * Xc1c2' / (Xc1c2*Xc1c2' + alpha*eye(nweights,nweights));
        yc1_train = Wc1 * Xc1c2;
        yc2_train = Wc2 * Xc1c2;
        votes(img_idxs, c1) = votes(img_idxs, c1) + (yc1_train>yc2_train)';
        votes(img_idxs, c2) = votes(img_idxs, c2) + (yc2_train>yc1_train)';
        W(:,c1,c2,1) = Wc1;
        W(:,c1,c2,2) = Wc2;
    end
end

[max_votes, classes] = max(votes');
nerrs = sum(test_classes(classes)~=Ytrain) 
err = nerrs / ntrain

% Test classifier
Ytest = Y(test_idxs);
Xtest = XX(widx, test_idxs);
votes = zeros(ntest, nclasses);
for c1=1:nclasses-1
    for c2=c1+1:nclasses
        yc1_test = W(:,c1,c2,1)' * Xtest;
        yc2_test = W(:,c1,c2,2)' * Xtest;
        votes(:, c1) = votes(:, c1) + (yc1_test>yc2_test)';
        votes(:, c2) = votes(:, c2) + (yc2_test>yc1_test)';
    end
end
[max_votes, classes] = max(votes');
nerrs = sum(test_classes(classes)~=Ytest) 
err = nerrs / ntest
