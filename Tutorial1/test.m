t = linspace(0,10, 100)';
y = sqrt(t).*linspace(1,1.5,5);
idx = find(t > 7, 1);
text_x = repmat(t(idx), 1, size(y,2));
text_y = y(idx, :);
labels = {'label1', 'label2', 'label3', 'label4', 'label5'};
plot(t,y)
hold on
t = text(text_x, text_y, labels, 'Color', 'k');