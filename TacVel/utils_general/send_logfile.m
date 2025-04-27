function send_logfile(emailto,attachedfilename)
%sends the "attachedfilename", e.g. the log file after running matlab
%analyses, to the "emailto" address, e.g. "vandortc@mit.edu".
%Chances are that your antivirus/firewall prevent this function from contacting SMPT server. 
%As a  workaound, you can temporarily turn off your antivirus/firewall. A
%better workaround would be to UNCHECK the "Scan outbound mail (SMTP)" 
%option in your antivirus settings. Or make MATLAB an exemption in your firewall settings.
%If using this function in the midst of a set of time taking calculation, toss it in a 
%'try catch' statement to prevent abortion of your program.

%   hozan@mit.edu 02/11/2016
%   settings for CJV lab matlabreport account, made exlusively for this
%   function. Ironically, yahoo was chosen over gmail because yahoo has often
%   less strict rules regarding SSL connections/security, or at least used
%   to be less strict.
mail        = 'matlabreport_cjvlab@yahoo.com';  %DO NOT use this email address for sensitive information as the password is written out loud in the next line.
password    = 'mouseemail';     
smtpserver  = 'smtp.mail.yahoo.com';
subject     = 'Matlab work is done!';
message     = sprintf(['\n\nLog file is attached:\t', attachedfilename,'\n' ,datestr(now)],'\n mohsen\n');

setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server',smtpserver);
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%% Send the email
sendmail(emailto,subject,message,attachedfilename)