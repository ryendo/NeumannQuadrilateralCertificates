classdef QuadrilateralProofRunner < handle
    % Unified entry point for the local and global certified computations.
    properties (SetAccess=private)
        Root
        IntlabRoot
        VeigsRoot
    end
    methods
        function self=QuadrilateralProofRunner(intlabRoot,veigsRoot)
            self.Root=fileparts(mfilename('fullpath'));
            if nargin<1, intlabRoot=''; end
            if nargin<2, veigsRoot=''; end
            self.IntlabRoot=intlabRoot;
            self.VeigsRoot=veigsRoot;
        end
        function setup(self)
            if ~isempty(self.IntlabRoot), addpath(self.IntlabRoot); end
            if ~isempty(self.VeigsRoot), addpath(self.VeigsRoot); end
            addpath(fullfile(self.Root,'src'));
            addpath(fullfile(self.Root,'tests'));
            if exist('startintlab','file')==2, startintlab; end
            assert(exist('intval','class')==8 || exist('intval','file')==2, ...
                'INTLAB is not available. Pass its root to QuadrilateralProofRunner.');
            assert(exist('veigs','file')==2 && exist('veig','file')==2, ...
                'veigs is not available. Pass its root to QuadrilateralProofRunner.');
        end
        function summary=summarizeLocal(self,resultsDir)
            if nargin<2, resultsDir=fullfile(self.Root,'results','local'); end
            summary=qn_summarize_local_results(resultsDir);
        end
        function runSingleBoxCertificateWorker(self,workerId,workerCount,faceSubdivisions,outDir)
            if nargin<5, outDir=fullfile(self.Root,'results','local_new'); end
            qn_local_certificate_cover(workerId,workerCount,faceSubdivisions,outDir);
        end
        function result=runGlobalCertifiedCover(self,outputFile,workerId,workerCount)
            if nargin<2, outputFile=fullfile(self.Root,'results','global','summary.json'); end
            if nargin<3, workerId=1; end
            if nargin<4, workerCount=1; end
            result=qn_global_certified_cover(3,60,true,outputFile,workerId,workerCount);
        end
        function report=smokeTest(self)
            report=qn_smoke_test(self.Root);
        end
    end
end
