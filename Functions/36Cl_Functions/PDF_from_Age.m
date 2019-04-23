function [ XPDFCosmo,PDFCosmo ] = PDF_from_Age( AgeCosmo,ErrCosmo )
%PDF_FROM_AGE Summary of this function goes here
%   Detailed explanation goes here

            % PDf boundaries
                NbSigma=7;
                AgeDiscr=0.01;% Determine time step

                AgeNSigmaInfBrut=AgeCosmo-NbSigma*ErrCosmo;
                AgeNSigmaSupBrut=AgeCosmo+NbSigma*ErrCosmo;
                
                if(AgeNSigmaInfBrut<=0); AgeNSigmaInfBrut=.0; end

                % Set on the grid
                AgeNSigmaInfDiscr=TrouvValeurDiscret(AgeNSigmaInfBrut,AgeDiscr);
                AgeNSigmaSupDiscr=TrouvValeurDiscret(AgeNSigmaSupBrut,AgeDiscr);

                % Create time vector for PDF
                XPDFCosmo=AgeNSigmaInfDiscr:AgeDiscr:AgeNSigmaSupDiscr;

                % Creat PDF
                
                PDFCosmo=normpdf(XPDFCosmo,AgeCosmo,ErrCosmo);
             
end

