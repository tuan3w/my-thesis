function plot(obj)

figure,
playbackPanelH = uipanel(gcf,'Units','Normalized','Position',[.3 0 .4 .1]);

buttWidth = 1/6;
centers = linspace(0,1,7)-buttWidth/2;
rewButtonH = uicontrol(playbackPanelH,'Style','pushbutton','String','<<','Units',...
	'Normalized','Position',[centers(2) 0.2 buttWidth 0.6]);
ffButtonH = uicontrol(playbackPanelH,'Style','togglebutton','String','>>','Units',...
	'Normalized','Position',[centers(3) 0.2 buttWidth 0.6]);
stopButtonH = uicontrol(playbackPanelH,'Style','pushbutton','String','stop','Units',...
	'Normalized','Position',[centers(4) 0.2 buttWidth 0.6]);
playButtonH = uicontrol(playbackPanelH,'Style','togglebutton','String','play','Units',...
	'Normalized','Position',[centers(5) 0.2 buttWidth 0.6]);
pauseButtonH = uicontrol(playbackPanelH,'Style','togglebutton','String','||','Units',...
	'Normalized','Position',[centers(6) 0.2 buttWidth 0.6]);

%% Plot the time domain signal
s = obj.s;
fs = obj.fs;
plot((1:length(s))/fs,s);
title('Audio signal')
xlabel('time (s)');
axis tight

player = audioplayer(s,fs);
set(player,'TimerPeriod',0.1);
set(player,'StartFcn',@plotTransportBar);
set(player,'TimerFcn',@updateTransportBar);
set(player,'StopFcn',@deleteTransportBar);

%% Add playback controls
set(playButtonH,'Callback',@play_callback);
set(stopButtonH,'Callback',@stop_callback);
set(pauseButtonH,'Callback',@pause_callback);
set(rewButtonH,'Callback',@rew_callback);
set(ffButtonH,'Callback',@ff_callback);

	function play_callback(~,~)
		set(player,'SampleRate',fs);
		play(player,player.CurrentSample);
		set(pauseButtonH,'Value',0);
		set(ffButtonH,'Value',0);
	end

	function pause_callback(~,~)
		pause(player);
		set(playButtonH,'Value',0);
		set(ffButtonH,'Value',0);
	end

	function stop_callback(~,~)
		stop(player);
		set(playButtonH,'Value',0);
		set(pauseButtonH,'Value',0);
		set(ffButtonH,'Value',0);
	end

	function ff_callback(~,~)
		set(player,'SampleRate',1.5*fs);
		set(pauseButtonH,'Value',0);
		set(playButtonH,'Value',0);
	end

	function rew_callback(~,~)
		stop(player);
		play(player);
		set(pauseButtonH,'Value',0);
		set(playButtonH,'Value',1);
	end

%% Transport Bar functions
	function plotTransportBar(~,~)
		global tbH
		xLim = get(gca,'Xlim');
		yLim = get(gca,'YLim');
		tbH = line([xLim(1) xLim(1)],yLim,'Color','k');
	end

	function updateTransportBar(hObject,~)
		global tbH
		currentSample = hObject.CurrentSample;
		pos = currentSample/fs;
		set(tbH,'XData',[pos pos]);
	end

	function deleteTransportBar(~,~)
		global tbH
		delete(tbH);
	end
end
