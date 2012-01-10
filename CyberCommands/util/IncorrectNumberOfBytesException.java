package util;

/*
 * Created on Mar 3, 2005
 *
 */


/**
 * @author wpun
 *
 */
public class IncorrectNumberOfBytesException extends Exception {

	private static final long serialVersionUID = 1L;

	public IncorrectNumberOfBytesException() {
		super();
	}

	public IncorrectNumberOfBytesException(String message) {
		super(message);
	}

	public IncorrectNumberOfBytesException(String message, Throwable cause) {
		super(message, cause);
	}

	public IncorrectNumberOfBytesException(Throwable cause) {
		super(cause);
	}

}
