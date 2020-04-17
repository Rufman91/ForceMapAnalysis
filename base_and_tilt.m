function [basedcurve,f,gof,deletecurve] = base_and_tilt(TH,vDef)
% subtract baseline and tilt from the forcecurve by fitting a function to
% the non contact domain the function tries to fit a non-affine-linear
% function. If the linear fit is too bad, the function tries a 9th grade
% polynomial fit instead

deletecurve = 0;
[ncd , ncidx] = no_contact_domain(vDef);
[f1,gof1] = fit(TH(1:ncidx),ncd,'poly1');

if gof1.rsquare < 0.9
    ftt = fittype( @(a,b,c,d,e,x) a+b*x+c*sin(d*x+e) );
    counter = 0;                                    
    gof2.rsquare = -1;
    Options = fitoptions('Method','NearestInterpolant');
    while gof2.rsquare < 0.9                            
        [f2 , gof2] = fit(TH(1:ncidx),ncd,ftt,Options);     
        counter = counter +1;
        if counter > 10
            f = f1;
            gof = gof1;
            warning("Wanted to but couldn't determine a sine-like-fit; taking linear fit instead")
            break
        else
        end
    end
    f = f2;
    gof = gof2;
else
    f = f1;
    gof = gof1;
end

if gof.rsquare < 0.9
        plot(TH,feval(f,TH),TH,vDef);
    answ = questdlg('Quality of base line fit is very low. consider excluding this force curve or retry fitting',...
        'Base line and tilt','Exclude it','Retry fit',...
        'Keep current fit','Keep current fit');
    if isequal(answ,'Exclude it')
        deletecurve = 1;
    elseif isequal(answ,'Retry fit')
        [basedcurve,f,gof,deletecurve] = base_and_tilt(vDef,TH);
        return
    else
    end
    close
end
basedcurve(:) = (vDef(:)-feval(f,TH));
end