int pwmpin = 3;


void setup() {
  // put your setup code here, to run once:
  pinMode(pwmpin, OUTPUT);
  analogWrite(pwmpin, 0);
  Serial.begin(9600);
  
}

void loop() {
  char cmd;
  char var;
  int val;
  String s;
  // Wait until someone prints something
  if (Serial.available() > 0)
  {
    cmd = Serial.read();
    if (cmd == '.'){
      delay(50);
      var = Serial.read();
      if (var == -1){
        Serial.println("ERROR - Command expected!");  
        delay(50);
        s = Serial.readString();
        return;
      }
      // Read the value
      delay(50);
      s = Serial.readStringUntil('\n');
      val = s.toInt();

      if (val > 255){
        val = 255;
      } else if (val < 0){
        val = 0;
      }
      analogWrite(pwmpin, val);
      
    }else{
      Serial.print("ERROR - Unknown command mode -->");
      Serial.println(cmd);
      delay(50);
      s = Serial.readString();
      return;
    }
  }

}
