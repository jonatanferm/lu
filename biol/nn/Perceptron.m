function gens = Perceptron
%Exercise written 9 August 2006 by Anders Brodin
%A one-cell perceptron with a bias and two inputs

%Global parameters:
global NUM_EXAMPLES NUM_INPUTS;                    
NUM_EXAMPLES = 4;
NUM_INPUTS = 2;
%Global arrays:
global INPUT CORRECT_OUTPUT W TOTAL_ERROR;
INPUT = zeros(NUM_EXAMPLES, NUM_INPUTS);
CORRECT_OUTPUT = zeros(1,NUM_EXAMPLES);
TOTAL_ERROR = zeros(1,NUM_EXAMPLES);
W = zeros(1,NUM_INPUTS + 1);                %W holds weights from inputs, NUM_IMPUTS + 1 is for bias weight

Initiate_Input;                             %function that initiates vector with input signals
Initiate_Expected;                          %initiates vector with correct output signals
Initiate_Weights;                           %initiates weights randomly between -5.0 and 5.0
eta = 0.8;                                  %coefficient for learning rate, close to 0 = low learning rate, close to 1 = high learning rate
i = 1;                                      %i must start as a positive number
Count = 0;                                  %counter for the number of batches (a batch = all four examples)
NUMTIMES = 0;
while i > 0 && NUMTIMES < 100                 %task will be solved when i = 0    
    for Example = 1 : NUM_EXAMPLES              %four examples of pairwise input 0 or 1
        SumIP = 0;                              %sums the product from the three weights * signals
        for Weight = 1 : NUM_INPUTS + 1         %Weight is counter for weights (2 inputs + bias)
            if Weight == NUM_INPUTS + 1         %if weight is bias
                Signal = 1;                     %get bias signal
            else
                Signal = INPUT(Example, Weight);%get signal from signal vector
            end
            IP = Signal * W(1, Weight);         %inner product
            SumIP = SumIP + IP;                 %sum for all weights
        end
        Output = StepFunction(SumIP);           %check if incoming signal produces output in transfer function      
        Error = Output - CORRECT_OUTPUT(1, Example);    %difference between expected and actual output (=supervision!)
        TOTAL_ERROR(1,Example) = abs(Error);            %vector where error for each example can be stored
        if Output > CORRECT_OUTPUT(1, Example)          %check if weights should be decreased
            W(1,1) = W(1,1) - abs(Error * INPUT(Example,1) * eta);      %adjust weight 1
            W(1,2) = W(1,2) - abs(Error * INPUT(Example,2) * eta);      %adjust weight 2
            W(1,3) = W(1,3) - abs(Error * eta);                         %adjust weight 3 (=bias)
        else
            if Output < CORRECT_OUTPUT(1, Example)      %check if weights should be increased
               W(1,1) = W(1,1) + abs(Error * INPUT(Example,1) * eta);   %adjust weight 1
               W(1,2) = W(1,2) + abs(Error * INPUT(Example,2) * eta);   %adjust weight 2
               W(1,3) = W(1,3) + abs(Error * eta);                      %adjust weight 3 (=bias)               
            else
                if (Example == NUM_EXAMPLES && sum(TOTAL_ERROR) == 0)   %if we are at example 4 with no error -> task is done, exit
                    i = 0;                                              %exit condition
                 end
            end
        end            
    end
    NUMTIMES = NUMTIMES + 1;
    Count = Count + 1;
    SaveOutput(Count) = sum(TOTAL_ERROR) / NUM_EXAMPLES;                %save mean error for plot          
    if i == 0 || NUMTIMES == 25
        %plot(SaveOutput)
        gens = length(SaveOutput);
    end
 end
%___________________________________________________
function Outsignal = StepFunction(Input)
%transfer function, decides whether summed input should produce output
if Input <= 0
    Outsignal = 0;
else
    Outsignal = 1;
end
%___________________________________________________
function Initiate_Weights
%initiates random numbers between -3 and 3 to weights
%input weighst are called W, output weight is called V
global W NUM_INPUTS;
for Count = 1 : NUM_INPUTS + 1
   if rand(1) < 0.5 
      W(Count) = 3 * rand(1);            %assign positive random number to W, vector for input weights
      if Count == 1
          V(1) = 3 * rand(1);            %positive ouput weight
      end
   else
      W(Count) = -3 * rand(1);           %assign negative random number to W, vector for input weights
      if Count == 1
          V(1) = -3 * rand(1);           %negative ouput weight
      end
   end
end

%___________________________________________________
function Initiate_Expected
%fills the vector CORRECT_OUTPUT with correct output...
%I start with expected output for "????"...
global CORRECT_OUTPUT;
CORRECT_OUTPUT(1,1) = 0;
CORRECT_OUTPUT(1,2) = 1;
CORRECT_OUTPUT(1,3) = 1;
CORRECT_OUTPUT(1,4) = 0;
%___________________________________________________
function Initiate_Input
%creates a vector "INPUT" for the four examples x two input neurons
global INPUT;
INPUT(1,1) = 0;             %input set 1
INPUT(1,2) = 0;
INPUT(2,1) = 1;             %inputset 2
INPUT(2,2) = 0;
INPUT(3,1) = 0;             %inputset 3
INPUT(3,2) = 1;
INPUT(4,1) = 1;             %inputset 4
INPUT(4,2) = 1;
%____________________________________________________
