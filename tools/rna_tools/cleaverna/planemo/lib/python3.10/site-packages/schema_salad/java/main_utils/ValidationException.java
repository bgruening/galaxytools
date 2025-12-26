package ${package}.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ValidationException extends RuntimeException {
  private final List<ValidationException> children;
  private String bullet = "";
  private String currentMessage;

  public ValidationException(final String message) {
    this(message, (List<ValidationException>) null);
  }

  public ValidationException(final String message, final ValidationException child) {
    this(message, Arrays.asList(child));
  }

  public ValidationException(final String message, final List<ValidationException> children_) {
    super(message);
    this.currentMessage = message;
    final List<ValidationException> children = new ArrayList<ValidationException>();
    if (children_ != null) {
      for (final ValidationException child : children_) {
        children.addAll(child.simplify());
      }
    }
    this.children = children;
  }

  public ValidationException withBullet(final String bullet) {
    this.bullet = bullet;
    return this;
  }

  public List<ValidationException> simplify() {
    if (getMessage().length() > 0) {
      return Arrays.asList(this);
    } else {
      return this.children;
    }
  }

  public String summary(final int level, final boolean withBullet) {
    final int indentPerLevel = 2;
    final String spaces = new String(new char[level * indentPerLevel]).replace("\0", " ");
    final String bullet;
    if (this.bullet.length() > 0 && withBullet) {
      bullet = this.bullet;
    } else {
      bullet = "";
    }
    return spaces + bullet + this.currentMessage;
  }

  public String prettyStr(final Integer level_) {
    Integer level = level_;
    if (level == null) {
      level = 0;
    }
    final List<String> parts = new ArrayList<String>();
    int nextLevel;
    if (this.currentMessage != null && this.currentMessage.length() > 0) {
      parts.add(this.summary(level, true));
      nextLevel = level + 1;
    } else {
      nextLevel = level;
    }
    for (final ValidationException child : this.children) {
      parts.add(child.prettyStr(nextLevel));
    }
    final String ret = String.join("\n", parts);
    return ret;
  }

  public String getMessage() {
    return this.prettyStr(null);
  }
}
