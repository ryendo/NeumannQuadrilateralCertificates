classdef QuadrilateralProofRunner < handle
    % Unified entry point for the local and global certified computations.
    properties (SetAccess=private)
        Root
        IntlabRoot
    end
    methods
        function self=QuadrilateralProofRunner(intlabRoot)
            self.Root=fileparts(mfilename('fullpath'));
            if nargin<1, intlabRoot=''; end
            self.IntlabRoot=intlabRoot;
        end
        function setup(self)
            if ~isempty(self.IntlabRoot), addpath(genpath(self.IntlabRoot)); end
            addpath(fullfile(self.Root,'src','local'));
            addpath(fullfile(self.Root,'src','global'));
            addpath(fullfile(self.Root,'tests'));
            if exist('startintlab','file')==2, startintlab; end
            assert(exist('intval','class')==8 || exist('intval','file')==2, ...
                'INTLAB is not available. Pass its root to QuadrilateralProofRunner.');
        end
        function summary=summarizeLocal(self,resultsDir)
            if nargin<2, resultsDir=fullfile(self.Root,'results','local'); end
            summary=qn_summarize_local_results(resultsDir);
        end
        function runLocalWorker(self,workerId,workerCount,faceSubdivisions,outDir)
            if nargin<5, outDir=fullfile(self.Root,'results','local_new'); end
            dq2_run_certificate(workerId,workerCount,faceSubdivisions,0,0,0,outDir);
        end
        function result=runGlobal(self,outputFile,workerId,workerCount)
            if nargin<2, outputFile=fullfile(self.Root,'results','global','summary.json'); end
            if nargin<3, workerId=1; end
            if nargin<4, workerCount=1; end
            result=qn_run_global_cover(3,60,true,outputFile,workerId,workerCount);
        end
        function report=smokeTest(self)
            report=qn_smoke_test(self.Root);
        end
    end
end
