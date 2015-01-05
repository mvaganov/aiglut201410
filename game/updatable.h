#pragma once

class Updatable {
public:
	virtual void update(const int ms) = 0;
};

class Transient {
public:
	virtual bool isAlive() const = 0;
	virtual bool isReadyToDelete() const = 0;
	virtual void markForDelete() = 0;
};