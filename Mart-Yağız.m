%% 3. Random Variable Y definition
close all;
clear all;
rand('seed', sum(100*clock));
num_samples = 1000000;

x1 = 2*rand(1,1000000);
x2 = -10*rand(1,1000000);
y_seq = x1+x2;

y_sqr = y_seq.*y_seq;
mean_sqr_samples = mean(y_sqr);
disp(mean_sqr_samples);

figure(1);
histogram(y_seq, 'Normalization', 'pdf');
xlabel('Y');
ylabel('Probability Density Function');
title('Histogram of Y');
grid on;

%% 4. Uniform Quantization 
N = 4;
a = zeros(1,N-1);
delta = 12/(N-1);
a(1)= -10 + delta/2;
for i = 1:(N-1-1)
    a(i+1) = a(i) + delta;
end

y_bar=zeros(1, num_samples);

for i = 1:length(y_seq)
    for j = 1:(N-1)
        if y_seq(i) <= a(j)
            y_bar(i) = -10 + delta*(j-1);
            break;
        elseif y_seq(i) > a(N-1)
            y_bar(i) = a(N-1) + delta/2;
        end
    end
end
q_err = y_seq - y_bar;
figure(2);
histogram(q_err);
xlabel('Quantization Error');
ylabel('Numbers');
title('Histogram of Quantization Error');
grid on;

q_error_sqr = q_err.*q_err;
avg_power = mean(q_error_sqr);
disp(avg_power);
SQNR_uniform = 10*log(mean_sqr_samples/avg_power);
disp(SQNR_uniform);

%% N=32 Quantization Levels
N = 32;
a = zeros(1,N-1);
delta = 12/(N-1);
a(1)= -10 + delta/2;
for i = 1:(N-1-1)
    a(i+1) = a(i) + delta;
end

y_bar = zeros(1, num_samples);

for i = 1:length(y_seq)
    for j = 1:(N-1)
        if y_seq(i) <= a(j)
            y_bar(i) = -10 + delta*(j-1);
            break;
        elseif y_seq(i) > a(N-1)
            y_bar(i) = a(N-1) + delta/2;
        end
    end
end
q_err = y_seq - y_bar;
figure(3);
histogram(q_err);
xlabel('Quantization Error');
ylabel('Numbers');
title('Histogram of Quantization Error');
grid on;

q_error_sqr = q_err.*q_err;
avg_power = mean(q_error_sqr);
disp(avg_power);
SQNR_uniform = 10*log(mean_sqr_samples/avg_power);
disp(SQNR_uniform);

%% Non-uniform Quantization:
syms fy(y);
fy(y)= @(y) piecewise((0<y<=2),(1/10 - y/20)*y,(-8<y<=0), y/10, (-8>=y>=-10), y*(y/20 + 1/2), (y<-10), 0, (y>2), 0);
syms f(y);
f(y)= @(y) piecewise((0<y<=2),1/10 - y/20,(-8<y<=0), 1/10, (-8>=y>=-10), y/20 + 1/2, (y<-10), 0, (y>2), 0);
disp(f(-7));

%% Lloyd max algorithm
%Step1:initial guess for reconstruction levels
N = 4;
a = zeros(1,N+1);
delta = 12/(N-1);
a(2)= -10 + delta/2;
for i = 2:(N-1)
    a(i+1) = a(i) + delta;
end
%% Step2: Finding optimal reconstruction levels
recon_levels = zeros(1,N);
a(1) = -Inf;
a(N+1) = Inf;
avg_power = [5.1,5];
y_bar2 = zeros(1, num_samples);
while (abs(avg_power(1) - avg_power(2))) > 0.00001
    for i = 1:N
        recon_levels(i) = int(fy(y),a(i),a(i+1))/int(f(y),a(i),a(i+1));
    end

    for h = 1:N-1
        a(h+1)=(recon_levels(h)+ recon_levels(h+1))/2;
    end
        for m = 1:length(y_seq)
            for j = 1:(N-1)
                if y_seq(m) <= a(j+1)
                    y_bar2(m) = recon_levels(j);
                    break;
                elseif y_seq(m) > a(N)
                    y_bar2(m) = recon_levels(N);
                end
            end
        end
    q_err = y_seq - y_bar2;
    q_error_sqr = q_err.*q_err;
    avg_power(1) = avg_power(2);
    avg_power(2) = mean(q_error_sqr);
end
reconstruction_indices = 1:length(recon_levels);

% Create a plot with markers for the reconstruction levels and vertical lines for boundaries
figure(8);
plot(reconstruction_indices, recon_levels, 'ro', 'MarkerSize', 8);
hold on;

% Plot horizontal lines connecting the reconstruction levels
for i = 1:length(recon_levels)
    if i < length(recon_levels)
        line([i, i + 1], [recon_levels(i), recon_levels(i)], 'Color', 'b');
    end
end

% Plot vertical lines at the boundaries
for i = 1:length(a(2:4))
    line([i+1, i+1], [recon_levels(i+1), recon_levels(i)], 'Color', 'b');
end

hold off;

% Customize the plot
title('Reconstruction Levels and Boundaries');
xlabel('Region Index');
ylabel('Reconstruction Levels');
xticks(1:length(a(2:4)));
xticklabels(cellstr(num2str(a(2:4)')));

grid on;
legend('Reconstruction Levels', 'Boundaries', 'Location', 'Best');


disp(avg_power(2));
SQNR_nonuniform = 10*log(mean_sqr_samples/avg_power(2));
disp(SQNR_nonuniform);
figure(4);
histogram(q_err);
xlabel('Quantization Error');
ylabel('Numbers');
title('Histogram of Quantization Error');
grid on;



%% N = 32 Non-uniform quantization
N = 32;
a = zeros(1,N+1);
delta = 12/(N-1);
a(2)= -10 + delta/2;
for i = 2:(N-1)
    a(i+1) = a(i) + delta;
end
recon_levels = zeros(1,N);
a(1) = -Inf;
a(N+1) = Inf;
avg_power = [5.1,5];
y_bar2 = zeros(1, num_samples);
count = 0;
while (abs(avg_power(1) - avg_power(2))) > 0.0001
%for p = 1:10
    count = count + 1;
    for i = 1:N
        recon_levels(i) = int(fy(y),a(i),a(i+1))/int(f(y),a(i),a(i+1));
    end

    for h = 1:N-1
        a(h+1)=(recon_levels(h)+ recon_levels(h+1))/2;
    end
        for m = 1:length(y_seq)
            for j = 1:(N-1)
                if y_seq(m) <= a(j+1)
                    y_bar2(m) = recon_levels(j);
                    break;
                elseif y_seq(m) > a(N)
                    y_bar2(m) = recon_levels(N);
                end
            end
        end
    q_err = y_seq - y_bar2;
    q_error_sqr = q_err.*q_err;
    avg_power(1) = avg_power(2);
    avg_power(2) = mean(q_error_sqr);
end
disp(count);
disp(avg_power(2));
SQNR_uniform = 10*log(mean_sqr_samples/avg_power(2));
disp(SQNR_uniform);
figure(5);
histogram(q_err);
xlabel('Quantization Error');
ylabel('Numbers');
title('Histogram of Quantization Error');
grid on;


%% 1. Lempel Ziv Encoding

text = 'At the dawn of time, the god known as the Creator created the Wheel of Time, which weaves the patterns of the times with the universe and the lives of men and women. The wheel has his seven spokes representing the era, and is turned by the only power that flows from the true source. One Power is divided into male and female halves, Saidin and Saidal, and they work together to drive the wheels while opposing each other. Humans who can harness that power are known as Channelers. Their main organization is called the Aes Sedai (an archaic term for all servants). The Creator imprisoned its antithesis, the dark one "Shaitan", at the moment of creation and kept it away from the wheel. However, in an era known as the Age of Legends, Aes Sedai experiments accidentally breached the Dark One prison, allowing his influence to once again permeate the world. He rallied the proud, the corrupt, and the ambitious to his cause, and these servants began efforts to free the dark Ones from their prisons once and for all, so that in his own image he could It made it possible to recreate reality. In response to this threat, the Wheel spun out dragons, channelers of immense power, to become champions of the Light. In Age of Legends, the dragon is a man named Lewes Serin Telamon, who eventually rose to command Ace Sedai and his allies in their battle against the forces of the Dark Ones. After ten years of grueling war, Lewes Serin leads his forces to victory in a daring attack on the Shayol Guul volcano (the site of the earthly link to the prison of the Dark One) and reseals the prison. was made. However, in a moment of triumph, the Dark Ones were able to defile Cydin, driving One Power male channelers insane. Lewes Seryn murdered her friends and family and committed suicide by deliberately overdosing on One Power. Other male channelers destroyed the world with One Power, causing earthquakes and tidal waves that reshaped the world. Eventually, the last male channeler was either killed or cut off from One Power, and humanity was all but extinct, leaving only women to safely wield One Power. Aes Sedai has rebuilt humanity and rescued it from this dark age. Mankind is now freed from prison by the Dark One, Dragons are reborn to fight him again, and even though he is Mankinds only hope against the Dark One, he destroys the world in an instant. lived under the shadow of the prophecy that it would. time in the process. Over the course of 3,000 years, humans have rebuilt technology to a level roughly equivalent to that of the late Middle Ages or early modern times. However, despite a higher level of general education and an understanding of hygiene and anatomy, there is a complete lack of formal scientific, industrial manufacturing, and academic institutions. This can be traced back to the nature of the post-apocalyptic world, where much knowledge survived but the structures and institutions that made that knowledge possible were destroyed.';
dict = {'#', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' '};
% Input text processing 
text = upper(text);
text = regexprep(text, '[^A-Z ]', ''); % Removes characters that are not uppercase letters or spaces
text = [text '#'];
% Initialize the output 
output = {};
% Encode Lempel-Ziv
i = 1;
prefix = text(i);

while i <= length(text)
    prefix = text(i);
    while i < length(text)
        if any(strcmp(dict, [prefix text(i+1)]))
            i = i + 1;
            prefix = append(prefix, text(i));
        else
            break;
        end
    end
% Append the binary representation of the index of the prefix to the output
    prefixIndex = find(strcmp(dict, prefix));
    output = [output, dec2bin(prefixIndex-1, ceil(log2(length(dict))))]; % -1 because MATLAB indices start at 1
% Add new string
    if i < length(text)
        i = i + 1;
        dict = [dict, append(prefix,text(i))];
    else

        break;
    end
end
%% 2. Calculate Average Wordlength
total_bit = 0;
bit_num = 0;

for i = 1:length(output)
    bit_num = cellfun('length',output(i));
    total_bit = total_bit + bit_num;
end
avg_codelength = total_bit / length(text);
disp(avg_codelength);


%% 3. Decode the output
dict2 = {'#', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' '};
decoded_text = '';
input_binary = output;
decimal_code = (bin2dec(input_binary) + 1)'; % Convert binary to decimal + 1 (to adjust for MATLAB indexing)
new_char = '';

previous_char = dict2(decimal_code(1));
for i=1:length(input_binary)
    new_char = '';
    message = decimal_code(i);
    if message <= length(dict2)
        new_char = dict2(message);
        decoded_text = append(decoded_text, new_char);
        if i>=2
        chr = cell2mat(new_char);
        added = append(previous_char,chr(1));
        dict2 = [dict2, added];
        previous_char = new_char;
        end
    else
        
    end
end

disp(cell2mat(decoded_text));

%% 4. Probability calculation
char_numbers = zeros(1,27);
dict3 = {'#', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' '};
for j = 2:length(dict3)
    i = 1;
    while i <= length(text)
        if strcmp(text(i),dict3(j))
            char_numbers(j-1) = char_numbers(j-1) + 1; 
            i = i+1;
        else
            i = i+1;
        end
    end
end
outcomes = 2:28;

% Calculate the total number of observations
total_observations = sum(char_numbers);

% Calculate the PMF by dividing frequencies by the total number of observations
pmf = char_numbers / total_observations;

% Create a bar plot for the PMF
bar(outcomes, pmf)

xlabel('Characters')
ylabel('Probability')
title('Probability Mass Function (PMF)')
xticks(1:28); % Assuming you have 26 data points (A to Z)
xticklabels({'#','A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' '});
grid on
pmf(11)=0;
avg_codeword = sum(pmf .* 5);
entropy = -sum(pmf .* log2(pmf + (pmf == 0)));

disp(entropy);





