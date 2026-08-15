// CNC Shield V3 - X & Z Move Only One Direction (120 mm)

#include <AccelStepper.h>

#define X_STEP_PIN 2
#define X_DIR_PIN  5
#define Z_STEP_PIN 4
#define Z_DIR_PIN  7
#define EN_PIN     8

AccelStepper stepperX(AccelStepper::DRIVER, X_STEP_PIN, X_DIR_PIN);
AccelStepper stepperZ(AccelStepper::DRIVER, Z_STEP_PIN, Z_DIR_PIN);

const long STEPS_TO_MOVE = 22000;   // 90 mm moving distance
const float MAX_SPEED = 1050.0;     // steps/sec
const float ACCEL = 200.0;          // steps/sec²

void setup()
{
  pinMode(EN_PIN, OUTPUT);
  digitalWrite(EN_PIN, LOW);   // Enable drivers

  stepperX.setMaxSpeed(MAX_SPEED);
  stepperX.setAcceleration(ACCEL);

  stepperZ.setMaxSpeed(MAX_SPEED);
  stepperZ.setAcceleration(ACCEL);

  // Move only once in one direction
  stepperX.move(STEPS_TO_MOVE);
  stepperZ.move(STEPS_TO_MOVE);
}

void loop()
{
  // Keep running until both motors reach the target
  if (stepperX.distanceToGo() != 0 || stepperZ.distanceToGo() != 0)
  {
    stepperX.run();
    stepperZ.run();
  }
  else
  {
    // Motion complete - stop here
    while (1);
  }
}