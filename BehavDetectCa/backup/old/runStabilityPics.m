%make average of recording images

strPath1 = 'D:\Data\Processed\imagingdata\20140715\xyt01\images';
intNumIm = round(10*25.4);


imTemp = zeros(512,512,3);
for intT=1:intNumIm
	for intCh=0:1
		strFile = sprintf('%s%st%05d_ch%02d.tif',strPath1,filesep,intT,intCh);
		imTemp(:,:,intCh+1)=imTemp(:,:,intCh+1)+im2double(imread(strFile));
	end
end
im1 = mean(cat(4,imTemp,circshift(imTemp,[1 0 0])),4);
vecG = buildGaussian(-1:1,0,1,0);
matG = vecG'*vecG;
im1(:,:,1) = conv2(im1(:,:,1),matG,'same');
im1(:,:,2) = conv2(im1(:,:,2),matG,'same');
im1 = imnorm(im1);
figure
imshow(im1)


strPath2 = 'D:\Data\Processed\imagingdata\20140715\xyt04\images';
intNumIm = round(10*25.4);

imTemp = zeros(512,512,3);
for intT=16400:(16400+intNumIm-1)
	for intCh=0:1
		strFile = sprintf('%s%st%05d_ch%02d.tif',strPath1,filesep,intT,intCh);
		imTemp(:,:,intCh+1)=imTemp(:,:,intCh+1)+im2double(imread(strFile));
	end
end
im2 = mean(cat(4,imTemp,circshift(imTemp,[1 0 0])),4);
vecG = buildGaussian(-1:1,0,1,0);
matG = vecG'*vecG;
im2(:,:,1) = conv2(im2(:,:,1),matG,'same');
im2(:,:,2) = conv2(im2(:,:,2),matG,'same');
im2 = imnorm(im2);
figure
imshow(im2)


%{
Appendix A. Algorithms for calculating adjusted p-values
There are various algorithms one could use to implement the procedures, especially
the computationally intensive ones. We describe algorithms for the computationally
intensive procedures that would be appropriate when there are large number of variables
but we are interested in only those with small adjusted p-values. For all the algorithms,
let the p-values of the observed data be (p1; p2; : : : ; pk ), and let the ordered p-values
for the observed data be p1 6p2 6· · ·6pk . The sequence {j} is 9xed throughout
the algorithms. Denote the adjusted p-values by (pˆ1 ;pˆ2 ; : : : ;pˆk ).
Algorithm A. This algorithm determines the adjusted p-values for Procedure A to allow
u false discoveries. The adjusted p-values are calculated recursively, starting by
setting pˆ1 = ˆ p2 = · · ·= ˆ pu = 0. The algorithm is used to calculate pˆr , the adjusted
p-value for variable r (r¿u) after having found pˆ1 ;pˆ2 ; : : : ;pˆr?1 .
Step 1. Choose B random permutations of the data vectors consistent with the experimental
design. Denote the univariate p-values for the variables from the jth permutation
by (p?( j)
1 ; p?( j)
2 ; : : : ; p?( j)
k ) for j = 1; 2; : : : ; B.
Step 2. Initialize COUNT =1. Let W ={w1; w2; : : : ; wu
} be a collection of distinct
indices from {1; 2; : : : ; r?1}. De9ne the set R = {w1; w2; : : : ; wu
} ? {r; r+1; : : : ; k}.
For each j = 1; 2; : : : ; B, 9nd "j , the (u + 1)st smallest p-value for the permuted data
over the subset R, i.e.,
"j = [{p?( j)
i : i ?R}](u+1);
where the notation [A](s) is de9ned to refer to the sth smallest of the elements of the
9nite set A.
Step 3. For j = 1; 2; : : : ; B, if pr ¿"j then COUNT ? COUNT + 1.
Step 4. p?
r (W) ? COUNT=(B + 1).
Step 5. For each W repeat steps 2–4. There will be ( r?1
u ) such sets W. (Note that
when u = 0, W can only be empty, and in that case R = {r; r+1; : : : ; k}.)
Step 6. p?
r
? maxW p?
r (W).
Step 7. pˆr
? max{pˆ1 ;pˆ2 ; : : : ;pˆr?1 ; ˆ p?r
}.
If the sample sizes are small enough so that all permutations of the data vectors
consistent with the experimental design can be enumerated, then j =1; : : : ; B in Step 1
should run over all permutations except the one corresponding to the observed data.
Algorithm B. This algorithm determines the adjusted p-values for Procedure B to allow
the false discovery proportion 6. The adjusted p-values are calculated recursively,
starting with pˆ1 . The algorithm below is used to calculate pˆr , the adjusted p-value for
variable r after having found pˆ1 ;pˆ2 ; : : : ;pˆr?1 . The algorithm is identical to Algorithm
A except rather than u being 9xed in advance, there is a Step 0:
Step 0. Let u=|[r]|. If u¿|[(r ?1)]|, then an automatic rejection is allowed and
set p?
r = 0 and skip to Step 7. Otherwise, continue with Step 1.
Algorithm A?. This algorithm determines the adjusted p-values for the conservative
version of Procedure A to allow u false discoveries. First, one starts by setting pˆ1 =
pˆ2=· · ·= ˆ pu=0. The algorithm then calculates all of the remaining adjusted p-values.
Step 0. Initialize the counters COUNTi = 1 for i = u + 1; : : : ; k.
Step 1. Choose a random permutation of the data vectors consistent with the experimental
design. Denote the univariate p-values for the variables from this permutation
by (p?
1; p?
2; : : : ; p?
k ).
Step 2. Let q? = [{p?
i : i = 1; : : : ; k}](u+1).
Step 3. If pi ¿q? then COUNTi ? COUNTi + 1 for i = u + 1; : : : ; k.
Step 4. Repeat Steps 1–3 B times.
Step 5. pˆi
? COUNTi=(B + 1) for i = u + 1; : : : ; k.
As with Algorithm A, if the sample sizes are small enough so that all permutations
of the data vectors consistent with the experimental design can be enumerated, then
all permutations except the one corresponding to the observed data should be used in
Steps 1–3.
Algorithm B?. This algorithm determines the adjusted p-values for the conservative
version of Procedure B to allow the false discovery proportion 6.
Step 0. Initialize the counters COUNTr = 1 for r = 1; : : : ; k.
Step 1. If |[r]|¿|[(r ? 1)]|, then p˜r = 0 for r = 1; : : : ; k:
Step 2. Choose a random permutation of the data vectors consistent with the experimental
design. Denote the univariate p-values for the variables from this permutation
by (p?
1; p?
2; : : : ; p?
k ). (Note that the ordering does not matter for this step.)
Step 3. If |[r]|=|[(r?1)]| and pr ¿[{p?
i : i=1; : : : ; k}](|[r]|+1) then COUNTr ?
COUNTr + 1 for r = 1; : : : ; k.
Step 4. Repeat Steps 2–3 B times.
Step 5. If |[r]| = |[(r ? 1)]| then p˜r
? COUNTr=(B + 1) for r = 1; : : : ; k.
Step 6. pˆr
? max(p˜1 ; : : : ;p˜r ) for r = 1; : : : ; k.
	%}