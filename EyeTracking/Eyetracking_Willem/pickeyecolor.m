function eyekleurpick_next = pickeyecolor(f,ex,x0,y0,video_heel,f1)
	
	g1 = round(f-(0.5*ex));
	g2 = round(f+(0.5*ex)-1);
	
	if g1<1
		
		g1=f;
		g2=f+ex-1;
		
	elseif f==f1
		
		g1=f;
		g2=f+10;
	end
	
	if g2>length(video_heel.frames)
		
		g1=length(video_heel.frames)-ex+1;
		g2=length(video_heel.frames);
		
	else
	end
	
	% for g=f:f+ex-1
	for g=g1:g2
		pick1 = im2double(video_heel.frames(1,g).cdata);
		[r c a] =size(pick1);
		img2 = zeros(r,c);
		pick = imadjust(pick1(:,:,1));
		%pick = pick1;
		
		    %Average pixels with surrounding, new matrix called img2
		
		    for m = 2:r-1
		        for n = 2:c-1
		            img2(m,n)=((pick(m,n)+pick(m,n+1)+pick(m,n-1)+pick(m+1,n)+pick(m+1,n+1)+pick(m+1,n-1)+pick(m-1,n)+pick(m-1,n+1)+pick(m-1,n-1))/9);
		        end
		    end
		
		    pick = img2;
		
% 		pi1(g) = pick(round(x0),round(y0));
% 		pi2(g) = pick(round(x0)+1,round(y0)+1);
% 		pi3(g) = pick(round(x0)-1,round(y0)-1);
% 		pi4(g) = pick(round(x0),round(y0)+1);
% 		pi5(g) = pick(round(x0),round(y0)-1);
% 		pi6(g) = pick(round(x0)+1,round(y0));
% 		pi7(g) = pick(round(x0)-1,round(y0));
% 		pi8(g) = pick(round(x0)-1,round(y0)+1);
% 		pi9(g) = pick(round(x0)+1,round(y0)-1);
% 		pi10(g) = pick(round(x0)+2,round(y0)+2);
% 		pi11(g) = pick(round(x0)-2,round(y0)-2);
% 		pi12(g) = pick(round(x0),round(y0)+2);
% 		pi13(g) = pick(round(x0),round(y0)-2);
% 		pi14(g) = pick(round(x0)+2,round(y0));
% 		pi15(g) = pick(round(x0)-2,round(y0));
% 		pi16(g) = pick(round(x0)-2,round(y0)+2);
% 		pi17(g) = pick(round(x0)+2,round(y0)-2);
% 		pi18(g) = pick(round(x0)-1,round(y0)+2);
% 		pi19(g) = pick(round(x0)-1,round(y0)-2);
% 		pi20(g) = pick(round(x0)-2,round(y0)-1);
% 		pi21(g) = pick(round(x0)-2,round(y0)+1);
% 		pi22(g) = pick(round(x0)+1,round(y0)+2);
% 		pi23(g) = pick(round(x0)+1,round(y0)-2);
% 		pi24(g) = pick(round(x0)+2,round(y0)-1);
% 		pi25(g) = pick(round(x0)+2,round(y0)+1);
		
		pi1(g) = pick(round(y0),round(x0));
		pi2(g) = pick(round(y0)+1,round(x0)+1);
		pi3(g) = pick(round(y0)-1,round(x0)-1);
		pi4(g) = pick(round(y0),round(x0)+1);
		pi5(g) = pick(round(y0),round(x0)-1);
		pi6(g) = pick(round(y0)+1,round(x0));
		pi7(g) = pick(round(y0)-1,round(x0));
		pi8(g) = pick(round(y0)-1,round(x0)+1);
		pi9(g) = pick(round(y0)+1,round(x0)-1);
		pi10(g) = pick(round(y0)+2,round(x0)+2);
		pi11(g) = pick(round(y0)-2,round(x0)-2);
		pi12(g) = pick(round(y0),round(x0)+2);
		pi13(g) = pick(round(y0),round(x0)-2);
		pi14(g) = pick(round(y0)+2,round(x0));
		pi15(g) = pick(round(y0)-2,round(x0));
		pi16(g) = pick(round(y0)-2,round(x0)+2);
		pi17(g) = pick(round(y0)+2,round(x0)-2);
		pi18(g) = pick(round(y0)-1,round(x0)+2);
		pi19(g) = pick(round(y0)-1,round(x0)-2);
		pi20(g) = pick(round(y0)-2,round(x0)-1);
		pi21(g) = pick(round(y0)-2,round(x0)+1);
		pi22(g) = pick(round(y0)+1,round(x0)+2);
		pi23(g) = pick(round(y0)+1,round(x0)-2);
		pi24(g) = pick(round(y0)+2,round(x0)-1);
		pi25(g) = pick(round(y0)+2,round(x0)+1);
		
	end
	
	pi1(pi1==0)=[];
	[f1a, f1b] = hist(pi1,20);f1a((f1b < 0.080)) = 0;
	ff(1) = f1b(find(f1a == max(f1a),1));
	pi2(pi2==0)=[];
	[f2a, f2b] = hist(pi2,20);f2a((f2b < 0.080)) = 0;
	ff(2) = f2b(find(f2a == max(f2a),1));
	pi3(pi3==0)=[];
	[f3a, f3b] = hist(pi3,20);f3a((f3b < 0.080)) = 0;
	ff(3) = f3b(find(f3a == max(f3a),1));
	pi4(pi4==0)=[];
	[f4a, f4b] = hist(pi4,20);f4a((f4b < 0.080)) = 0;
	ff(4) = f4b(find(f4a == max(f4a),1));
	pi5(pi5==0)=[];
	[f5a, f5b] = hist(pi5,20);f5a((f5b < 0.080)) = 0;
	ff(5) = f5b(find(f5a == max(f5a),1));
	pi6(pi6==0)=[];
	[f6a, f6b] = hist(pi6,20);f6a((f6b < 0.080)) = 0;
	ff(6) = f6b(find(f6a == max(f6a),1));
	pi7(pi7==0)=[];
	[f7a, f7b] = hist(pi7,20);f7a((f7b < 0.080)) = 0;
	ff(7) = f7b(find(f7a == max(f7a),1));
	pi8(pi8==0)=[];
	[f8a, f8b] = hist(pi8,20);f8a((f8b < 0.080)) = 0;
	ff(8) = f8b(find(f8a == max(f8a),1));
	pi9(pi9==0)=[];
	[f9a, f9b] = hist(pi9,20);f9a((f9b < 0.080)) = 0;
	ff(9) = f9b(find(f9a == max(f9a),1));
	pi10(pi10==0)=[];
	[f10a, f10b] = hist(pi10,20); f10a((f10b < 0.080)) = 0;
	ff(10) = f10b(find(f10a == max(f10a),1));
	pi11(pi11==0)=[];
	[f11a, f11b] = hist(pi11,20); f11a((f11b < 0.080)) = 0;
	ff(11) = f11b(find(f11a == max(f11a),1));
	pi12(pi12==0)=[];
	[f12a, f12b] = hist(pi12,20); f12a((f12b < 0.080)) = 0;
	ff(12) = f12b(find(f12a == max(f12a),1));
	pi13(pi13==0)=[];
	[f13a, f13b] = hist(pi13,20); f13a((f13b < 0.080)) = 0;
	ff(13) = f13b(find(f13a == max(f13a),1));
	pi14(pi14==0)=[];
	[f14a, f14b] = hist(pi14,20); f14a((f14b < 0.080)) = 0;
	ff(14) = f14b(find(f14a == max(f14a),1));
	pi15(pi15==0)=[];
	[f15a, f15b] = hist(pi15,20); f15a((f15b < 0.080)) = 0;
	ff(15) = f15b(find(f15a == max(f15a),1));
	pi16(pi16==0)=[];
	[f16a, f16b] = hist(pi16,20); f16a((f16b < 0.080)) = 0;
	ff(16) = f16b(find(f16a == max(f16a),1));
	pi17(pi17==0)=[];
	[f17a, f17b] = hist(pi17,20); f17a((f17b < 0.080)) = 0;
	ff(17) = f17b(find(f17a == max(f17a),1));
	pi18(pi18==0)=[];
	[f18a, f18b] = hist(pi18,20); f18a((f18b < 0.080)) = 0;
	ff(18) = f18b(find(f18a == max(f18a),1));
	pi19(pi19==0)=[];
	[f19a, f19b] = hist(pi19,20); f19a((f19b < 0.080)) = 0;
	ff(19) = f19b(find(f19a == max(f19a),1));
	pi20(pi20==0)=[];
	[f20a, f20b] = hist(pi20,20); f20a((f20b < 0.080)) = 0;
	ff(20) = f20b(find(f20a == max(f20a),1));
	pi21(pi21==0)=[];
	[f21a, f21b] = hist(pi21,20); f21a((f21b < 0.080)) = 0;
	ff(21) = f21b(find(f21a == max(f21a),1));
	pi22(pi22==0)=[];
	[f22a, f22b] = hist(pi22,20); f22a((f22b < 0.080)) = 0;
	ff(22) = f22b(find(f22a == max(f22a),1));
	pi23(pi23==0)=[];
	[f23a, f23b] = hist(pi23,20); f23a((f23b < 0.080)) = 0;
	ff(23) = f23b(find(f23a == max(f23a),1));
	pi24(pi24==0)=[];
	[f24a, f24b] = hist(pi24,20); f24a((f24b < 0.080)) = 0;
	ff(24) = f24b(find(f24a == max(f24a),1));
	pi25(pi25==0)=[];
	[f25a, f25b] = hist(pi25,20); f25a((f25b < 0.080)) = 0;
	ff(25) = f25b(find(f25a == max(f25a),1));
	
	%maken van vector van alle eyekleurpicks, dan een histogram waarbij
	%outliers vervangen worden door gemiddelde over bepaalde frame.
	
	ff(ff<0.1)=[];
	sort(ff,'descend');
	if length(ff)>2
		ff=ff(1:3);
	end
	eyekleurpick_next = abs(mean(ff));
