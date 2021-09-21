# Train a model.
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

iris = load_iris()
X, y = iris.data, iris.target
X, y = X[:20], y[:20]
loss = list()
X_train, X_test, y_train, y_test = train_test_split(X, y)
clr = RandomForestClassifier(n_estimators=5)
clr.fit(X_train, y_train)
clr.score(X_test, y_test)
